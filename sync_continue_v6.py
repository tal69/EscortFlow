#!/usr/bin/env python3
"""Pull missing/newer continue_v6 files from bukhan1 using the Linux access skill."""

import argparse
import base64
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile


ALIAS = "bukhan1"
DESTINATION = Path(__file__).resolve().parent / "results_continue_v6"
HELPER = Path.home() / ".codex/skills/lab-linux-access/scripts/lab_access.py"
REMOTE_MANIFEST = r'''python3 - <<'PY'
import json
from pathlib import Path
root = (Path.home() / "EscortFlow/results_continue_v6").resolve(strict=True)
files = {}
for path in sorted(root.rglob("*")):
    if path.is_symlink():
        raise RuntimeError("Unexpected remote symlink: " + str(path))
    if path.is_file():
        stat = path.stat()
        files[path.relative_to(root).as_posix()] = {
            "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "mode": stat.st_mode & 0o777,
        }
print(json.dumps({"root": str(root), "files": files}))
PY'''
REMOTE_SNAPSHOTS = r'''
import base64
import json
from pathlib import Path
import sys
root = Path.home() / "EscortFlow/results_continue_v6"
snapshots = {}
for name in json.load(sys.stdin):
    path = root / name
    for attempt in range(5):
        before = path.stat()
        data = path.read_bytes()
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns):
            snapshots[name] = {
                "data": base64.b64encode(data).decode("ascii"),
                "size": len(data), "mtime_ns": after.st_mtime_ns,
                "mode": after.st_mode & 0o777,
            }
            break
    else:
        raise RuntimeError("Could not read a consistent snapshot: " + name)
print(json.dumps(snapshots))
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--campaign", choices=("continue_v6", "continue_v6_six_targets"),
                        default="continue_v6", help="original campaign or the revised six-target campaign")
    args = parser.parse_args()
    campaign_directory = "results_" + args.campaign
    destination_root = DESTINATION.with_name(campaign_directory)
    remote_manifest = REMOTE_MANIFEST.replace("EscortFlow/results_continue_v6", "EscortFlow/" + campaign_directory)
    remote_snapshots = REMOTE_SNAPSHOTS.replace("EscortFlow/results_continue_v6", "EscortFlow/" + campaign_directory)
    spec = importlib.util.spec_from_file_location("lab_access", HELPER)
    lab = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lab)
    machine = lab.machine(ALIAS)
    env = lab.authentication_env(machine)
    ssh = [lab.SSH, *lab.ssh_options(machine), "-T", "-l", machine["user"],
           machine["ssh_host"]]
    result = subprocess.run([*ssh, remote_manifest], env=env,
                            stdin=subprocess.DEVNULL, capture_output=True,
                            text=True, check=True, timeout=60, start_new_session=True)
    if result.stderr:
        print(result.stderr, end="")
    manifest = json.loads(result.stdout)
    selected = {}
    missing = 0
    for name, metadata in manifest["files"].items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeError("Unsafe remote relative path: " + name)
        destination = destination_root / relative
        if any(path.is_symlink() for path in (destination, *destination.parents)):
            raise RuntimeError("Unexpected destination symlink: " + str(destination))
        if destination.exists() and not destination.is_file():
            raise RuntimeError("Destination is not a regular file: " + str(destination))
        if not destination.exists():
            missing += 1
            selected[name] = metadata
        elif metadata["mtime_ns"] > destination.stat().st_mtime_ns:
            selected[name] = metadata
    print(f"Source: {machine['ssh_host']}:{manifest['root']}")
    print(f"Destination: {destination_root}")
    print(f"Remote files: {len(manifest['files'])}; missing: {missing}; "
          f"newer: {len(selected) - missing}; skipped: {len(manifest['files']) - len(selected)}")
    if args.dry_run:
        for name in selected:
            print("Would copy: " + name)
        return
    if not selected:
        print("Already up to date.")
        return
    destination_root.mkdir(parents=True, exist_ok=True)
    # Stage on the same filesystem so replacement is atomic. Never delete local files.
    with tempfile.TemporaryDirectory(prefix=".sync_continue_v6-", dir=destination_root.parent) as temp:
        staging = Path(temp) / "files"
        staging.mkdir()
        file_list = Path(temp) / "files.list"
        file_list.write_bytes(b"".join(os.fsencode(name) + b"\0" for name in selected))
        remote_shell = shlex.join([lab.SSH, *lab.ssh_options(machine)])
        command = [shutil.which("rsync") or "/usr/bin/rsync", "-rt", "--from0",
                   "--files-from=" + str(file_list), "--stats", "-e", remote_shell,
                   lab.remote_path(machine, manifest["root"] + "/"), str(staging) + "/"]
        result = lab.execute(command, env, 300)
        if result:
            raise RuntimeError(f"rsync failed with exit code {result}; nothing was installed.")
        # Capture rapidly growing outputs in one remote read, with before/after
        # timestamps checked there instead of across multiple network round trips.
        changed = [name for name, metadata in selected.items()
                   if (staging / name).stat().st_size != metadata["size"]]
        if changed:
            result = subprocess.run(
                [*ssh, "python3 -c " + shlex.quote(remote_snapshots)], env=env,
                input=json.dumps(changed), capture_output=True, text=True,
                check=True, timeout=60, start_new_session=True)
            if result.stderr:
                print(result.stderr, end="")
            snapshots = json.loads(result.stdout)
            for name in changed:
                metadata = snapshots[name]
                data = base64.b64decode(metadata.pop("data"), validate=True)
                if len(data) != metadata["size"]:
                    raise RuntimeError("Snapshot size verification failed: " + name)
                (staging / name).write_bytes(data)
                selected[name] = metadata
                print("Captured active file: " + name)
        for name, metadata in selected.items():
            if (staging / name).stat().st_size != metadata["size"]:
                raise RuntimeError("Transferred size verification failed: " + name)
        copied = 0
        total_bytes = 0
        for name, metadata in selected.items():
            destination = destination_root / name
            if destination.exists() and destination.stat().st_mtime_ns >= metadata["mtime_ns"]:
                print("Skipped local file updated during transfer: " + name)
                continue
            source = staging / name
            source.chmod(metadata["mode"])
            # Protocol 29 truncates mtimes to seconds. Restore full source precision.
            os.utime(source, ns=(source.stat().st_atime_ns, metadata["mtime_ns"]))
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, destination)
            stat = destination.stat()
            if stat.st_size != metadata["size"] or stat.st_mtime_ns != metadata["mtime_ns"]:
                raise RuntimeError("Destination size/timestamp verification failed: " + name)
            copied += 1
            total_bytes += stat.st_size
            print("Copied: " + name)
        print(f"Verified {copied} copied files ({total_bytes:,} bytes).")


if __name__ == "__main__":
    main()
