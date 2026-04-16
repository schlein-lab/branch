# BRANCH Code Review: Classify Layer (Teil 2/4)

**Reviewer**: service-zyrkel  
**Datum**: 2025-01-21  
**Scope**: `src/classify/{features.hpp, cascade.hpp}`

## Files

| File | LOC | Zweck |
|------|-----|-------|
| features.hpp | 42 | 12 Feature-Enums, FeatureVector (48B float array), BubbleCandidate struct |
| cascade.hpp | 85 | CascadeConfig, 4-Stage classify_one() mit early-exit |

## Critical Issues

| Datei:Zeile | Severity | Problem |
|-------------|----------|----------|
| cascade.hpp:61 | MEDIUM | `stage_depth_ratio()` gibt raw DepthRatio zurück, aber Threshold ist 1.8 — Semantik unklar: ist 1.8 = 180% Coverage? Sollte normalisiert sein [0,1] |
| cascade.hpp:78-80 | HIGH | Stage 2 setzt confidence=1.0f hart statt p2 zu nutzen — Confidence ist nicht kalibriert |
| cascade.hpp:90-92 | MEDIUM | Stage 4 nutzt SegdupAnnotationFlag als Prior, aber ignoriert RepeatAnnotationFlag — asymmetrisch |
| cascade.hpp:37-39 | LOW | CascadeConfig::min_bubble_length_bp/min_coverage/max_recursion_depth werden NICHT in classify_one() genutzt — Dead Config |

## Classification Logic Feedback

**Positiv:**
- ✅ Early-exit Cascade-Pattern korrekt implementiert
- ✅ 4 Stages in sinnvoller Reihenfolge: Flanken → Depth → Span → Prior
- ✅ Feature-Enum stabil für ML-Input-Index
- ✅ BubbleCandidate hat genomische Koordinaten (chrom_id, start, end)

**Logik-Probleme:**
- ❌ **Stage 2 Threshold-Semantik**: DepthRatio ≥ 1.8 → Duplication ist biologisch falsch. Duplikation hat Coverage ~2x diploid, also DepthRatio ~2.0. Threshold 1.8 ist zu niedrig (false positives bei normalen Regionen mit 1.8x Coverage durch Mappability-Artefakte)
- ❌ **Keine Mixed-Klasse-Erkennung**: classify_one() gibt nie `BubbleClass::Mixed` zurück — der dokumentierte Recursion-Mechanismus fehlt komplett
- ❌ **Confidence nicht kalibriert**: Stages geben teils 1.0f, teils raw Feature-Werte als Confidence — nicht vergleichbar

**Fehlende Discriminatoren:**
- ❌ ReadSpanCoverageIQR (Feature 3) wird nicht genutzt
- ❌ RefAlignScoreLeft/Right (Features 4,5) werden nicht genutzt
- ❌ PathLengthDifference, GcContentDivergence (Features 6,7) werden nicht genutzt
- ❌ AlleleFrequencyEstimate (Feature 10) wird nicht genutzt
- ❌ BubbleLengthBp (Feature 11) wird nicht genutzt

## TODOs

1. **MIXED_CLASS**: Implementiere rekursive Dekomposition für Mixed-Fälle mit den 5 Termination-Garantien aus dem Header-Kommentar
2. **DEPTH_THRESHOLD**: Ändern auf ~2.0 (oder ML-trained in v0.3) mit SegDup-Ausnahme
3. **CALIBRATION**: Alle Stage-Confidences auf [0,1] normalisieren via Platt-Scaling o.ä.
4. **USE_ALL_FEATURES**: Features 3-7, 10-11 in Cascade einbauen (oder für LightGBM in v0.3 dokumentieren)
5. **DEAD_CONFIG**: min_bubble_length_bp etc. tatsächlich in classify_one() prüfen

---
*Review 2/4 abgeschlossen. Nächster Scope: backend/backend_vtable.hpp, gpu/overlap_kernel.cuh*
