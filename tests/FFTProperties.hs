{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeApplications #-}
-- This module uses the test-framework-quickcheck2 package.
module Main where

import Control.Monad
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Storable as VS
import Data.Complex

import Test.Framework (defaultMain, testGroup)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import qualified Numeric.FFT.Vector.Invertible as I
import qualified Numeric.FFT.Vector.Invertible.Multi as IM
import qualified Numeric.FFT.Vector.Unitary as U
import qualified Numeric.FFT.Vector.Unitary.Multi as UM
import Numeric.FFT.Vector.Plan

main = defaultMain
            -- NB: There's no explicit tests for the Unnormalized package.
            -- However, its Planners are implicitly used by the other modules,
            -- so it's covered in the below tests.
            [ testGroup "invertibility"
              [ testProperty "I.dft" $ prop_invert (I.dft @Double) I.idft
              , testProperty "I.dftR2C" $ prop_invert (I.dftR2C @Double) I.dftC2R
              , testProperty "I.dct1" $ prop_invert (I.dct1 @Double) I.idct1
              , testProperty "I.dct2" $ prop_invert (I.dct2 @Double) I.idct2
              , testProperty "I.dct3" $ prop_invert (I.dct3 @Double) I.idct3
              , testProperty "I.dct4" $ prop_invert (I.dct4 @Double) I.idct4
              , testProperty "I.dst1" $ prop_invert (I.dst1 @Double) I.idst1
              , testProperty "I.dst2" $ prop_invert (I.dst2 @Double) I.idst2
              , testProperty "I.dst3" $ prop_invert (I.dst3 @Double) I.idst3
              , testProperty "I.dst4" $ prop_invert (I.dst4 @Double) I.idst4
              , testProperty "U.dft" $ prop_invert (U.dft @Double) U.idft
              , testProperty "U.dftR2C" $ prop_invert (U.dftR2C @Double) U.dftC2R
              , testProperty "U.dct2" $ prop_invert (U.dct2 @Double) U.idct2
              ]
            , testGroup "orthogonality"
              [ testProperty "U.dft" $ prop_orthog (U.dft @Double)
              , testProperty "U.idft" $ prop_orthog (U.idft @Double)
              , testProperty "U.dftR2C" $ prop_orthog (U.dftR2C @Double)
              , testProperty "U.dftC2R" $ prop_orthog (U.dftR2C @Double)
              , testProperty "U.dct2" $ prop_orthog (U.dct2 @Double)
              , testProperty "U.idct2" $ prop_orthog (U.idct2 @Double)
              , testProperty "U.dct4" $ prop_orthog (U.dct4 @Double)
              ]
            , testGroup "invertibility ND"
              [ testProperty "IM.dft" $ prop_invertND (IM.dft @Double) IM.idft
              , testProperty "IM.dftR2C" $ prop_invertND (IM.dftR2C @Double) IM.dftC2R
              , testProperty "UM.dft" $ prop_invertND (UM.dft @Double) UM.idft
              , testProperty "UM.dftR2C" $ prop_invertND (UM.dftR2C @Double) UM.dftC2R
              ]
            , testGroup "orthogonality"
              [ testProperty "UM.dft" $ prop_orthogND (UM.dft @Double)
              , testProperty "UM.idft" $ prop_orthogND (UM.idft @Double)
              , testProperty "UM.dftR2C" $ prop_orthogND (UM.dftR2C @Double)
              , testProperty "UM.dftC2R" $ prop_orthogND (UM.dftR2C @Double)
              ]
            ]

-------------------
-- An instance of Arbitrary that probably belongs in another package.

instance (V.Unbox a, Arbitrary a) => Arbitrary (V.Vector a) where
    arbitrary = V.fromList `fmap` arbitrary


-------------------------
-- Support functions to compare Doubles for (near) equality.

class Num a => Mag a where
    mag :: a -> Double

instance Mag Double where
    mag = abs

instance Mag (Complex Double) where
    mag = magnitude

-- Robustly test whether two Doubles are nearly identical.
close :: Mag a => a -> a -> Bool
close x y = tol > mag (x-y) / max 1 (mag x + mag y)
  where
    tol = 1e-10

withinTol :: (Mag a, V.Unbox a) => V.Vector a -> V.Vector a -> Bool
withinTol a b
    | V.length a /= V.length b = False
    | otherwise = V.and $ V.zipWith close a b


---------------------
-- The actual properties

-- Test whether the inverse actually inverts the forward transform.
prop_invert f g a = let
                        p1 = plan f (V.length a)
                        p2 = plan g (V.length a)
                    in (V.length a > 1) ==> withinTol a $ execute p2 $ execute p1 a

-- Test whether the transform preserves the L2 (sum-of-squares) norm.
prop_orthog f a = let
                    p1 = plan f (V.length a)
                  in (V.length a > 1) ==> close (norm2 a) (norm2 $ execute p1 a)

data DimsAndValues a = DimsAndValues (VS.Vector Int) (V.Vector a)
  deriving (Show)

instance (Arbitrary a, V.Unbox a) => Arbitrary (DimsAndValues a) where
  arbitrary = do
    dims <- liftM (VS.fromList . map getPositive) arbitrary `suchThatMap` maybeReduceSize
    values <- V.replicateM (VS.product dims) arbitrary
    return (DimsAndValues dims values)
    where
      -- We use this to prevent test cases from growing too big
      maybeReduceSize ds =
        if VS.product ds < 1000 then Just ds else maybeReduceSize (VS.init ds)

prop_invertND f g (DimsAndValues ds a) = let
                        p1 = planND f ds
                        p2 = planND g ds
                    in (V.length a > 1) ==> withinTol a $ execute p2 $ execute p1 a

prop_orthogND f (DimsAndValues ds a) = let
                    p1 = planND f ds
                  in (V.length a > 1) ==> close (norm2 a) (norm2 $ execute p1 a)

norm2 a = sqrt $ V.sum $ V.map (\x -> x*x) $ V.map mag a
