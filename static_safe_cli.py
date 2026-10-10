"""Current benchmark options for the two standalone static solver commands.

This only translates option names and fixes the formulation. The shared runner
executes the selected solver; it does not launch a campaign or another process.
"""
ALIASES = {
    '--Lx': '-x', '--Ly': '-y', '--output_cells': '-O',
    '--escorts_range': '-e', '--load_num': '-l', '--reps_range': '-r',
    '--retrieval_mode': '-m', '--csv': '-f', '--num_threads': '--threads',
    '-t': '--weighted-time-limit', '--time_limit': '--weighted-time-limit',
    '--total_time_limit': '--weighted-time-limit',
}


def run(formulation, argv):
    from RunSafeWeightedStatic import main
    translated = []
    for argument in argv:
        option, separator, value = argument.partition('=')
        if option == '--safe-weighted':
            if separator:
                raise SystemExit('--safe-weighted takes no value')
            continue
        if option == '--formulation':
            raise SystemExit('The standalone command fixes the formulation; omit --formulation')
        # The current benchmark always uses Gurobi and the common greedy start.
        if option in ('--gurobi', '--warmstart') and not separator:
            continue
        translated.append(ALIASES.get(option, option) + (separator + value if separator else ''))
    return main(['--formulation', formulation, *translated])
