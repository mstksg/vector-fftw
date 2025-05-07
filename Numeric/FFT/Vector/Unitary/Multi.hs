{- |
This module provides normalized versions of the transforms in @fftw@.

All of the transforms are normalized so that

 - Each transform is unitary, i.e., preserves the inner product and the sum-of-squares norm of its input.

 - Each backwards transform is the inverse of the corresponding forwards transform.

(Both conditions only hold approximately, due to floating point precision.)

For more information on the underlying transforms, see
<http://www.fftw.org/fftw3_doc/What-FFTW-Really-Computes.html>.
--
-- @since 0.2
-}

module Numeric.FFT.Vector.Unitary.Multi
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

import Control.Exception (assert)
import Control.Monad (forM_)
import Numeric.FFT.Vector.Base
import Numeric.FFT.Vector.FFI
import qualified Numeric.FFT.Vector.Unnormalized.Multi as U
import Data.Complex
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as MS
import Control.Monad.Primitive(RealWorld)

-- | A discrete Fourier transform. The output and input sizes are the same (@n@).
--
-- @y_k = (1\/sqrt n) sum_(j=0)^(n-1) x_j e^(-2pi i j k\/n)@
dft :: FFTW a => TransformND a (Complex a) (Complex a)
dft = U.dft {normalizationND = \ns -> constMultOutput scaleComplex $ 1 / sqrt (toEnum (VS.product ns))}

-- | An inverse discrete Fourier transform.  The output and input sizes are the same (@n@).
--
-- @y_k = (1\/sqrt n) sum_(j=0)^(n-1) x_j e^(2pi i j k\/n)@
idft :: FFTW a => TransformND a (Complex a) (Complex a)
idft = U.idft {normalizationND = \ns -> constMultOutput scaleComplex $ 1 / sqrt (toEnum (VS.product ns))}

-- | A forward discrete Fourier transform with real data.  If the input size is @n@,
-- the output size will be @n \`div\` 2 + 1@.
dftR2C :: FFTW a => TransformND a a (Complex a)
dftR2C = base {normalizationND = \ns -> modifyOutput $
                    complexR2CScaling (sqrt 2) ns (outputSizeND base $ VS.last ns)
        }
  where
    base = U.dftR2C

-- | A normalized backward discrete Fourier transform which is the left inverse of
-- 'U.dftR2C'.  (Specifically, @run dftC2R . run dftR2C == id@.)
--
-- This 'Transform' behaves differently than the others:
--
--  - Calling @plan dftC2R n@ creates a 'Plan' whose /output/ size is @n@, and whose
--    /input/ size is @n \`div\` 2 + 1@.
--
--  - If @length v == n@, then @length (run dftC2R v) == 2*(n-1)@.
--
dftC2R :: FFTW a => TransformND a (Complex a) a
dftC2R = base {normalizationND = \ns -> modifyInput $
                    complexR2CScaling (sqrt 0.5) ns (inputSizeND base $ VS.last ns)
        }
  where
    base = U.dftC2R


complexR2CScaling :: (MS.Storable a, Floating a, Enum a) => a -> VS.Vector Int -> Int -> MS.MVector RealWorld (Complex a) -> IO ()
complexR2CScaling !t !ns !len !a = assert (MS.length a == VS.product (VS.init ns) * len) $ do
    let !s1 = sqrt (1/toEnum (VS.product ns))
    let !s2 = t * s1
    -- Justification for the use of unsafeModify:
    -- The output size is 2n+1; so if n>0 then the output size is >=1;
    -- and if n even then the output size is >=3.
    forM_ [0.. VS.product (VS.init ns) - 1] $ \idx -> do
      unsafeModify a (idx * len) $ scaleComplex s1
      if odd (VS.last ns)
        then multC scaleComplex s2 (MS.unsafeSlice (idx * len + 1) (len-1) a)
        else do
            unsafeModify a (idx * len + len - 1) $ scaleComplex s1
            multC scaleComplex s2 (MS.unsafeSlice (idx * len + 1) (len-2) a)

