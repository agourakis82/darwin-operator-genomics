# Results: Quaternionic Compression Hypothesis

## Experimental Setup

- **N chains**: 1000 (800 train, 200 test)
- **Chain length**: 50 operators
- **Seed**: 42

## Results

### Group Case (Invertible Operators Only: S, R, K, RC, M)

| Model | Accuracy | Perplexity |
|-------|----------|------------|
| Markov(1) | 20.6% | 5.0 |
| Markov(2) | 20.2% | 5.01 |
| Quaternion | 19.5% | 5.04 |

### Semigroup Case (With Indels: D, I, V)

| Model | Accuracy | Perplexity |
|-------|----------|------------|
| Markov(1) | 11.9% | 8.01 |
| Markov(2) | 12.6% | 8.05 |
| Quaternion | 12.3% | 8.02 |

## Interpretation

**Negative Result**: The quaternion model does not consistently outperform Markov baselines.

This suggests that simple quaternionic representations may not capture meaningful structure
in random operator chains. Possible explanations:

1. The quaternion model lacks sufficient capacity or training
2. Random operator chains lack the sequential structure that quaternions might capture
3. The quaternionic hypothesis may not apply to this domain

This negative result is preserved for scientific completeness. Future work could explore:
- More sophisticated quaternion architectures
- Real evolutionary operator chains (not random)
- Alternative algebraic representations (Clifford algebras, etc.)

---
*Generated: 2025-12-14T19:12:54.873*
