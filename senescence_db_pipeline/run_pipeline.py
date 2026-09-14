"""End-to-end pipeline runner. Executes steps 1-8 in order.

Usage:
    python run_pipeline.py            # run all steps
    python run_pipeline.py 4 5 6      # run only the listed steps (reuses prior outputs)
"""
import sys
from pipeline import (step1_acquire, step2_harmonize, step3_unify, step4_conservation,
                      step5_expression, step6_interactions, step7_assemble, step8_metrics)

STEPS = {
    1: step1_acquire, 2: step2_harmonize, 3: step3_unify, 4: step4_conservation,
    5: step5_expression, 6: step6_interactions, 7: step7_assemble, 8: step8_metrics,
}

def main():
    want = [int(x) for x in sys.argv[1:]] or list(STEPS)
    for i in want:
        print(f"\n===== STEP {i}: {STEPS[i].__name__} =====")
        STEPS[i].run()
    print("\n[done] pipeline complete. See results/ and figures/.")

if __name__ == "__main__":
    main()
