# Limitations

## Model limitations

- The results are based on a reduced-order model integrating 1.5D fluid-structure behavior and 2D truth-backed closure, not full 3D offshore validation.
- The critical velocity Ucrit is currently reconstructed from tau_w, particle radius, center concentration, and wall concentration because the closure table does not contain an explicit Ucrit column.
- Segregation amplitude saturates at the imposed upper bound, suggesting that the high-segregation closure should be revisited experimentally or numerically.
- Long-duration deposition management remains incomplete.

## Diagnostic limitations

- Residual pressure-gradient time-derivative robustness was tested only against additive Gaussian noise, one-sample delay, and 5-point smoothing.
- Sensor drift, missing data, bias, calibration errors, and harsh field-measurement conditions remain untested.

## Engineering limitations

- No sea-trial validation.
- No full offshore engineering design.
- No environmental permitting analysis.
- No industrial safety certification.
- No field-ready operational procedure.

## Repository limitations

- Virtual environment and compiled dependencies are excluded.
- This archive is a paper companion, not a certified reproduction environment or full deployment package.
