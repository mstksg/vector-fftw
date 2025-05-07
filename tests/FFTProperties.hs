{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications #-}
-- This module uses the test-framework-quickcheck2 package.
module Main where

import Control.Monad
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Unboxed as V
import Data.Complex
import Data.Proxy

import Test.Framework (defaultMain, testGroup, Test)
import Test.Framework.Providers.QuickCheck2 (testProperty)
import Test.QuickCheck

import qualified Numeric.FFT.Vector.Invertible as I
import qualified Numeric.FFT.Vector.Invertible.Multi as IM
import qualified Numeric.FFT.Vector.Unitary as U
import qualified Numeric.FFT.Vector.Unitary.Multi as UM
import Numeric.FFT.Vector.Plan

main = defaultMain
            [ testGroup "double" (suite (Proxy @Double))
            , testGroup "float" (suite (Proxy @Float))
            ]

-- NB: There's no explicit tests for the Unnormalized package.
-- However, its Planners are implicitly used by the other modules,
-- so it's covered in the below tests.
suite
    :: forall a. (FFTW a, VG.Vector V.Vector a, Mag a a, Mag (Complex a) a, V.Unbox a, Arbitrary a, Show a) => Proxy a -> [Test]
suite p =   [ testGroup "invertibility"
              [ testProperty "I.dft" $ prop_invert (I.dft @a) I.idft
              , testProperty "I.dftR2C" $ prop_invert (I.dftR2C @a) I.dftC2R
              , testProperty "I.dct1" $ prop_invert (I.dct1 @a) I.idct1
              , testProperty "I.dct2" $ prop_invert (I.dct2 @a) I.idct2
              , testProperty "I.dct3" $ prop_invert (I.dct3 @a) I.idct3
              , testProperty "I.dct4" $ prop_invert (I.dct4 @a) I.idct4
              , testProperty "I.dst1" $ prop_invert (I.dst1 @a) I.idst1
              , testProperty "I.dst2" $ prop_invert (I.dst2 @a) I.idst2
              , testProperty "I.dst3" $ prop_invert (I.dst3 @a) I.idst3
              , testProperty "I.dst4" $ prop_invert (I.dst4 @a) I.idst4
              , testProperty "U.dft" $ prop_invert (U.dft @a) U.idft
              , testProperty "U.dftR2C" $ prop_invert (U.dftR2C @a) U.dftC2R
              , testProperty "U.dct2" $ prop_invert (U.dct2 @a) U.idct2
              ]
            , testGroup "orthogonality"
              [ testProperty "U.dft" $ prop_orthog (U.dft @a)
              , testProperty "U.idft" $ prop_orthog (U.idft @a)
              , testProperty "U.dftR2C" $ prop_orthog (U.dftR2C @a)
              , testProperty "U.dftC2R" $ prop_orthog (U.dftR2C @a)
              , testProperty "U.dct2" $ prop_orthog (U.dct2 @a)
              , testProperty "U.idct2" $ prop_orthog (U.idct2 @a)
              , testProperty "U.dct4" $ prop_orthog (U.dct4 @a)
              ]
            , testGroup "invertibility ND"
              [ testProperty "IM.dft" $ prop_invertND (IM.dft @a) IM.idft
              , testProperty "IM.dftR2C" $ prop_invertND (IM.dftR2C @a) IM.dftC2R
              , testProperty "UM.dft" $ prop_invertND (UM.dft @a) UM.idft
              , testProperty "UM.dftR2C" $ prop_invertND (UM.dftR2C @a) UM.dftC2R
              ]
            , testGroup "orthogonality"
              [ testProperty "UM.dft" $ prop_orthogND (UM.dft @a)
              , testProperty "UM.idft" $ prop_orthogND (UM.idft @a)
              , testProperty "UM.dftR2C" $ prop_orthogND (UM.dftR2C @a)
              , testProperty "UM.dftC2R" $ prop_orthogND (UM.dftR2C @a)
              ]
            ]

-------------------
-- An instance of Arbitrary that probably belongs in another package.

instance (V.Unbox a, Arbitrary a) => Arbitrary (V.Vector a) where
    arbitrary = V.fromList `fmap` arbitrary


-------------------------
-- Support functions to compare Doubles for (near) equality.

class (Num a, Floating b, Ord b) => Mag a b | a -> b where
    mag :: a -> b
    tol :: b

instance Mag Double Double where
    mag = abs
    tol = 1e-10

instance Mag (Complex Double) Double where
    mag = magnitude
    tol = 1e-10

instance Mag Float Float where
    mag = abs
    tol = 1e-4

instance Mag (Complex Float) Float where
    mag = magnitude
    tol = 1e-4

-- Robustly test whether two Doubles are nearly identical.
close :: forall a b. Mag a b => a -> a -> Bool
close x y = tol @a > mag (x-y) / max 1 (max (mag x) (mag y))

withinTol :: (Mag a b, V.Unbox a) => V.Vector a -> V.Vector a -> Bool
withinTol a b
    | V.length a /= V.length b = False
    | otherwise = V.and $ V.zipWith close a b


---------------------
-- The actual properties

-- Test whether the inverse actually inverts the forward transform.
prop_invert f g a = let
                        p1 = plan f (V.length a)
                        p2 = plan g (V.length a)
                        a' = execute p2 $ execute p1 a
                    in counterexample (show a') $ (V.length a > 1) ==> withinTol a a'

-- Test whether the transform preserves the L2 (sum-of-squares) norm.
prop_orthog f a = let
                    p1 = plan f (V.length a)
                    b = execute p1 a
                  in counterexample (show (norm2 a, norm2 b)) $ (V.length a > 1) ==> close (norm2 a) (norm2 b)

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
