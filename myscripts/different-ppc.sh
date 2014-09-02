#!/bin/bash

(cd ../../../output; echo ppc=30; ../projects/particles/bin/peano-Particles pit video 250000 30 0.05 0.1)
(cd ../../../output; echo ppc=50; ../projects/particles/bin/peano-Particles pit video 250000 50 0.05 0.1)
(cd ../../../output; echo ppc=100; ../projects/particles/bin/peano-Particles pit video 250000 100 0.05 0.1)
