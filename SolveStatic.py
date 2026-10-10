#!/usr/bin/env python3
"""Solve static SBM instances with either formulation and the corrected protocol.

Select --formulation loadflow or escortflow. Both choices automatically use the
same greedy warm start, sufficient integer objective coefficient, and physical
horizon. See README.md for examples and RunContinue.sh for the matched campaign.
"""

from RunSafeWeightedStatic import main


if __name__ == "__main__":
    raise SystemExit(main())
