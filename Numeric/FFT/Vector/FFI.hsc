
-- | All FFI functions, to unify things under a single typeclass
module Numeric.FFT.Vector.FFI(
              FFTW(..),
              CPlan(..),
              CFlags,
              CDirection,
              CKind,
              fftForward,
              fftBackward,
              PlanType(..),
              Preservation(..),
              planInitFlags,
            ) where

import Data.Complex
import Data.Word

import Foreign (Storable(..), Ptr, FunPtr, ForeignPtr, (.|.))
import Foreign.C (CInt(..), CUInt(..))

#include <fftw3.h>

type CFlags = CUInt

newtype CPlan b = CPlan {unCPlan :: ForeignPtr (CPlan b) }

foreign import ccall unsafe "fftw_execute" fftw_execute_double :: Ptr (CPlan Double) -> IO ()
foreign import ccall "&fftw_destroy_plan" fftw_destroy_plan_double :: FunPtr (Ptr (CPlan Double) -> IO ())

foreign import ccall unsafe "fftwf_execute" fftw_execute_float :: Ptr (CPlan Float) -> IO ()
foreign import ccall "&fftwf_destroy_plan" fftw_destroy_plan_float :: FunPtr (Ptr (CPlan Float) -> IO ())

-- | Whether the complex fft is forwards or backwards.
type CDirection = CInt

fftForward :: CDirection
fftForward = (#const FFTW_FORWARD)

fftBackward :: CDirection
fftBackward = (#const FFTW_BACKWARD)

---------------------
-- Creating FFTW plans

-- First, the Transform flags:
data PlanType = Estimate | Measure | Patient | Exhaustive
data Preservation = PreserveInput | DestroyInput

-- | Marshal the Transform flags for use by fftw.
planInitFlags :: PlanType -> Preservation -> CFlags
planInitFlags pt pr = planTypeInt .|. preservationInt
  where
    planTypeInt = case pt of
                    Estimate -> #const FFTW_ESTIMATE
                    Measure -> #const FFTW_MEASURE
                    Patient -> #const FFTW_PATIENT
                    Exhaustive -> #const FFTW_EXHAUSTIVE
    preservationInt = case pr of
                    PreserveInput -> #const FFTW_PRESERVE_INPUT
                    DestroyInput -> #const FFTW_DESTROY_INPUT


-- | The type of the cosine or sine transform.
type CKind = (#type fftw_r2r_kind)

foreign import ccall unsafe "fftw_plan_dft_1d" fftw_plan_dft_1d_double
    :: CInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> CDirection
        -> CFlags -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_dft_1d" fftw_plan_dft_1d_float
    :: CInt -> Ptr (Complex Float) -> Ptr (Complex Float) -> CDirection
        -> CFlags -> IO (Ptr (CPlan Float))

foreign import ccall unsafe "fftw_plan_dft_r2c_1d" fftw_plan_dft_r2c_1d_double
    :: CInt -> Ptr Double -> Ptr (Complex Double) -> CFlags -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_dft_r2c_1d" fftw_plan_dft_r2c_1d_float
    :: CInt -> Ptr Float -> Ptr (Complex Float) -> CFlags -> IO (Ptr (CPlan Float))

foreign import ccall unsafe "fftw_plan_dft_c2r_1d" fftw_plan_dft_c2r_1d_double
    :: CInt -> Ptr (Complex Double) -> Ptr Double -> CFlags -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_dft_c2r_1d" fftw_plan_dft_c2r_1d_float
    :: CInt -> Ptr (Complex Float) -> Ptr Float -> CFlags -> IO (Ptr (CPlan Float))

foreign import ccall unsafe "fftw_plan_r2r_1d" fftw_plan_r2r_1d_double
    :: CInt -> Ptr Double -> Ptr Double -> CKind -> CFlags -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_r2r_1d" fftw_plan_r2r_1d_float
    :: CInt -> Ptr Float -> Ptr Float -> CKind -> CFlags -> IO (Ptr (CPlan Float))

foreign import ccall unsafe "fftw_plan_dft" fftw_plan_dft_double
    :: CInt -> Ptr CInt -> Ptr (Complex Double) -> Ptr (Complex Double)
        -> CDirection -> CFlags -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_dft" fftw_plan_dft_float
    :: CInt -> Ptr CInt -> Ptr (Complex Float) -> Ptr (Complex Float)
        -> CDirection -> CFlags -> IO (Ptr (CPlan Float))

foreign import ccall unsafe "fftw_plan_dft_r2c" fftw_plan_dft_r2c_double
    :: CInt -> Ptr CInt -> Ptr Double -> Ptr (Complex Double) -> CFlags
        -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_dft_r2c" fftw_plan_dft_r2c_float
    :: CInt -> Ptr CInt -> Ptr Float -> Ptr (Complex Float) -> CFlags
        -> IO (Ptr (CPlan Float))

foreign import ccall unsafe "fftw_plan_dft_c2r" fftw_plan_dft_c2r_double
    :: CInt -> Ptr CInt -> Ptr (Complex Double) -> Ptr Double -> CFlags
        -> IO (Ptr (CPlan Double))

foreign import ccall unsafe "fftwf_plan_dft_c2r" fftw_plan_dft_c2r_float
    :: CInt -> Ptr CInt -> Ptr (Complex Float) -> Ptr Float -> CFlags
        -> IO (Ptr (CPlan Float))

class (Storable a, Floating a, Enum a) => FFTW a where
    fftw_execute :: Ptr (CPlan a) -> IO ()
    fftw_destroy_plan :: FunPtr (Ptr (CPlan a) -> IO ())
    fftw_plan_dft_1d :: CInt -> Ptr (Complex a) -> Ptr (Complex a) -> CDirection -> CFlags -> IO (Ptr (CPlan a))
    fftw_plan_dft_r2c_1d :: CInt -> Ptr a -> Ptr (Complex a) -> CFlags -> IO (Ptr (CPlan a))
    fftw_plan_dft_c2r_1d :: CInt -> Ptr (Complex a) -> Ptr a -> CFlags -> IO (Ptr (CPlan a))
    fftw_plan_r2r_1d :: CInt -> Ptr a -> Ptr a -> CKind -> CFlags -> IO (Ptr (CPlan a))
    fftw_plan_dft :: CInt -> Ptr CInt -> Ptr (Complex a) -> Ptr (Complex a) -> CDirection -> CFlags -> IO (Ptr (CPlan a))
    fftw_plan_dft_r2c :: CInt -> Ptr CInt -> Ptr a -> Ptr (Complex a) -> CFlags -> IO (Ptr (CPlan a))
    fftw_plan_dft_c2r :: CInt -> Ptr CInt -> Ptr (Complex a) -> Ptr a -> CFlags -> IO (Ptr (CPlan a))

instance FFTW Double where
    fftw_execute = fftw_execute_double
    fftw_destroy_plan = fftw_destroy_plan_double
    fftw_plan_dft_1d = fftw_plan_dft_1d_double
    fftw_plan_dft_r2c_1d = fftw_plan_dft_r2c_1d_double
    fftw_plan_dft_c2r_1d = fftw_plan_dft_c2r_1d_double
    fftw_plan_r2r_1d = fftw_plan_r2r_1d_double
    fftw_plan_dft = fftw_plan_dft_double
    fftw_plan_dft_r2c = fftw_plan_dft_r2c_double
    fftw_plan_dft_c2r = fftw_plan_dft_c2r_double

instance FFTW Float where
    fftw_execute = fftw_execute_float
    fftw_destroy_plan = fftw_destroy_plan_float
    fftw_plan_dft_1d = fftw_plan_dft_1d_float
    fftw_plan_dft_r2c_1d = fftw_plan_dft_r2c_1d_float
    fftw_plan_dft_c2r_1d = fftw_plan_dft_c2r_1d_float
    fftw_plan_r2r_1d = fftw_plan_r2r_1d_float
    fftw_plan_dft = fftw_plan_dft_float
    fftw_plan_dft_r2c = fftw_plan_dft_r2c_float
    fftw_plan_dft_c2r = fftw_plan_dft_c2r_float

