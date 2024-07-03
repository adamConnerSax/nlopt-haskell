{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE CApiFFI #-}

{- |
Module      :  Numeric.Optimization.NLOPT.NNLS
Copyright   :  (c) Adam Conner-Sax
License     :  BSD3
Maintainer  :  Adam Conner-Sax
Stability   :  provisional
Portability :  GHC

Low-level interface to the NNSL, LDP and LSI functions inside the  NLOPT library.
-}

module Numeric.Optimization.NLOPT.NNLS where

import Numeric.Optimization.NLOPT.Bindings as NB

import Foreign hiding (void)
import Foreign.C.String
import Foreign.C.Types
import qualified Foreign.Concurrent as CFP

import qualified Data.Vector.Storable.Mutable as MV
import qualified Data.Vector.Storable as V

-- These routines expect matrices to be column major.
-- Hmatrix "flatten" concatenated rows.
-- So we flatten the transpose, thus storing it as column major and the leading dimension
-- is the number of rows.

foreign import ccall "nnls_"
  nnls_ :: Ptr CDouble -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CDouble
           -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble -> Ptr CInt
           -> Ptr CInt -> IO ()

-- argmin_x ||Ax - b||_2 s.t. x >= 0
nnls :: Int -> Int -> V.Vector Double -> V.Vector Double -> IO (V.Vector Double, Double, Int)
nnls m n aV bV = do
  let lda = m -- transposed row major
  aMV <- V.thaw $ V.unsafeCast aV
  bMV <- V.thaw $ V.unsafeCast bV
  MV.unsafeWith aMV $ \a ->
    with (fromIntegral lda) $ \ldaPtr ->
      with (fromIntegral n) $ \nPtr ->
       with (fromIntegral m) $ \mPtr ->
         MV.unsafeWith bMV $ \b ->
           allocaBytes (cdblBytes 1) $ \rNorm ->
             allocaBytes (cdblBytes n) $ \w ->
               allocaBytes (cdblBytes m) $ \z ->
                 allocaBytes (cintBytes n) $ \indx ->
                   allocaBytes (cintBytes 1) $ \mode -> do
                     xMV <- MV.unsafeNew n
                     MV.unsafeWith xMV $ \x -> do
                       nnls_ a ldaPtr mPtr nPtr b x rNorm w z indx mode
                     xV <- V.unsafeCast <$> V.freeze xMV
                     mode' <- fromIntegral <$> peek mode
                     rNorm' <- realToFrac <$> peek rNorm
                     pure (xV, rNorm', mode')

foreign import ccall "ldp_"
  ldp_ :: Ptr CDouble -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble -> Ptr CInt -> Ptr CInt -> IO ()


-- argmin_x ||x||_2 s.t. Gx >= h
ldp :: Int -> Int -> V.Vector Double -> V.Vector Double -> IO (V.Vector Double, Double, Int)
ldp m n gV hV = do
  let ldg = m
      g' = V.unsafeCast gV
      h' = V.unsafeCast hV
  V.unsafeWith g' $ \gPtr ->
    with (fromIntegral ldg) $ \ldgPtr ->
      with (fromIntegral m) $ \mPtr ->
        with (fromIntegral n) $ \nPtr ->
          V.unsafeWith h' $ \hPtr -> do
            xMV <- MV.unsafeNew n
            MV.unsafeWith xMV $ \xPtr ->
              allocaBytes (cdblBytes 1) $ \xNormPtr -> do
              let xSize = (m + 2) * (n + 1) + (2 * m)
              allocaBytes (cdblBytes xSize) $ \wPtr ->
                allocaBytes (cintBytes m) $ \indxPtr ->
                  allocaBytes (cintBytes 1) $ \modePtr -> do
                    ldp_ gPtr ldgPtr mPtr nPtr hPtr xPtr xNormPtr wPtr indxPtr modePtr
                    xV <- V.unsafeCast <$> V.freeze xMV
                    mode <- fromIntegral <$> peek modePtr
                    xNorm <- realToFrac <$> peek xNormPtr
                    pure (xV, xNorm, mode)



foreign import ccall "lsi_"
  lsi_ :: Ptr CDouble -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble
       -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CInt
       -> Ptr CDouble -> Ptr CDouble -> Ptr CDouble
       -> Ptr CInt -> Ptr CInt -> IO ()

-- x is length n, which is also # cols of E and G
lsi :: Int -> Int -> Int -> V.Vector Double -> V.Vector Double -> V.Vector Double -> V.Vector Double -> IO (V.Vector Double, Double, Int)
lsi n me mg e f g h = do
  let le = me
      lg = mg
  eMV <- V.thaw $ V.unsafeCast e
  fMV <- V.thaw $ V.unsafeCast f
  gMV <- V.thaw $ V.unsafeCast g
  hMV <- V.thaw $ V.unsafeCast h
  MV.unsafeWith eMV $ \ePtr ->
    MV.unsafeWith fMV $ \fPtr ->
      MV.unsafeWith gMV $ \gPtr ->
        MV.unsafeWith hMV $ \hPtr ->
          with (fromIntegral le) $ \lePtr ->
            with (fromIntegral me) $ \mePtr ->
              with (fromIntegral lg) $ \lgPtr ->
                with (fromIntegral mg) $ \mgPtr ->
                  with (fromIntegral n) $ \nPtr -> do
                    xMV <- MV.unsafeNew n
                    MV.unsafeWith xMV $ \xPtr ->
                      allocaBytes (cdblBytes 1) $ \xNormPtr -> do
                        let wSize = (n + 1) * (mg + 2) + 2 * mg
                        allocaBytes (cdblBytes wSize) $ \wPtr ->
                          allocaBytes (cintBytes lg) $ \jwPtr ->
                            allocaBytes (cintBytes 1) $ \modePtr -> do
                             lsi_ ePtr fPtr gPtr hPtr lePtr mePtr lgPtr mgPtr nPtr xPtr xNormPtr wPtr jwPtr modePtr
                             xV <- V.unsafeCast <$> V.freeze xMV
                             mode <- fromIntegral <$> peek modePtr
                             xNorm <- realToFrac <$> peek xNormPtr
                             pure (xV, xNorm, mode)

cdblBytes :: Int -> Int
cdblBytes k = sizeOf (undefined :: CDouble) * k

cintBytes :: Int -> Int
cintBytes k = sizeOf (undefined :: CInt) * k
--nnls_(aMV, aCols', aCols', aRows', bMV, xMV, rNorm', wMV, zMV', indxMV, mode')
