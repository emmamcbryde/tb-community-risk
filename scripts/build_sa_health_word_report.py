from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.apy.sa_health_report import build_sa_health_word_report


def main() -> None:
    result = build_sa_health_word_report()
    print("SA Health APY health-economic Word report generated.")
    print(f"Output directory: {result['outputDir']}")
    print(f"Word report: {result['docxPath']}")
    print(f"PDF verification: {result['pdfPath']}")
    print(f"Page count: {result['pageCount']}")


if __name__ == "__main__":
    main()
