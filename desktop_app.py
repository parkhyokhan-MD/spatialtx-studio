from __future__ import annotations

import sys

from spatialtx_desktop import __version__


def main() -> int:
    if any(argument in {"-h", "--help"} for argument in sys.argv[1:]):
        print(f"SpatialTX Studio Desktop v{__version__}")
        print("Usage: spatialtx-desktop")
        print("Launches the Windows desktop application. Main Mapper analyzes H5AD; raw inputs use Import / Convert.")
        return 0
    from spatialtx_desktop.app import main as desktop_main

    desktop_main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
