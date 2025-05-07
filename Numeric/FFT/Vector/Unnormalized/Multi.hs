{- |
Raw, unnormalized multi-dimensional versions of the transforms in @fftw@.

Note that the forwards and backwards transforms of this module are not actually
inverses.  For example, @run idft (run dft v) /= v@ in general.

For more information on the individual transforms, see
<http://www.fftw.org/fftw3_doc/What-FFTW-Really-Computes.html>.

@since 0.2
-}

module Numeric.FFT.Vector.Unnormalized.Multi
  (
        -- * Creating and executing 'Plan's
        run,
        plan,
        execute,
        -- * Complex-to-complex transforms
        dft,
        idft,
        -- * Real-to-complex transforms
        dftR2C,
        dftC2R,
  ) where

import Numeric.FFT.Vector.Base
import Numeric.FFT.Vector.FFI
import Data.Complex

dftND :: FFTW a => CDirection -> TransformND a (Complex a) (Complex a)
dftND d = TransformND
  { inputSizeND = id
  , outputSizeND = id
  , creationSizeFromInputND = id
  , makePlanND = \rk dims a b -> withPlanner . fftw_plan_dft rk dims a b d
  , normalizationND = const id
  }

-- | A forward discrete Fourier transform.  The output and input sizes are the same (@n@).
dft :: FFTW a => TransformND a (Complex a) (Complex a)
dft = dftND fftForward

-- | A backward discrete Fourier transform.  The output and input sizes are the same (@n@).
idft :: FFTW a => TransformND a (Complex a) (Complex a)
idft = dftND fftBackward

-- | A forward discrete Fourier transform with real data.  If the input size is @n0 * ... * nk@,
-- the output size will be @n0 * ... * nk \`div\` 2 + 1@.
dftR2C :: FFTW a => TransformND a a (Complex a)
dftR2C = TransformND {
              inputSizeND = id,
              outputSizeND = \n -> n `div` 2 + 1,
              creationSizeFromInputND = id,
              makePlanND = \rk dims a b -> withPlanner . fftw_plan_dft_r2c rk dims a b,
              normalizationND = const id
          }

-- | A backward discrete Fourier transform which produces real data.
--
-- This 'Transform' behaves differently than the others:
--
--  - Calling @plan dftC2R n@ creates a 'Plan' whose /output/ size is @n@, and whose
--    /input/ size is @n \`div\` 2 + 1@.
--
--  - If @length v == n@, then @length (run dftC2R v) == 2*(n-1)@.
dftC2R :: FFTW a => TransformND a (Complex a) a
dftC2R = TransformND {
            inputSizeND = \n -> n `div` 2 + 1,
            outputSizeND = id,
            creationSizeFromInputND = \n -> 2 * (n-1),
            makePlanND = \rk dims a b -> withPlanner . fftw_plan_dft_c2r rk dims a b,
            normalizationND = const id
        }
