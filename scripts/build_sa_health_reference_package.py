from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.apy.sa_health_reference_package import write_sa_health_reference_package


def main() -> None:
    result = write_sa_health_reference_package()
    manifest = result["manifest"]
    print("SA Health APY working-reference package generated.")
    print(f"Output directory: {result['outputDir']}")
    print(f"Configuration hash: {manifest['configurationHash']}")
    print(f"Economics configuration hash: {manifest['economicsConfigurationHash']}")
    print(f"Evidence registry hash: {manifest['evidenceRegistryHash']}")
    print("Interpretation: working-default analysis for planning - provisional.")


if __name__ == "__main__":
    main()
