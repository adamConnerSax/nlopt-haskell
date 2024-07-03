{-# OPTIONS_GHC -fno-warn-name-shadowing #-}
{-# OPTIONS_GHC -Wall #-}

module Main where

import Numeric.Optimization.NLOPT.Bindings
import Numeric.Optimization.NLOPT.NNLS
import qualified Data.Vector.Storable as V
import qualified System.Exit as E
import Control.Monad (when)

objective :: ScalarFunction ()
objective params g _ = do
  case g of
    Just grad -> V.copy grad $ V.fromList [0.0, 0.5]
    Nothing -> return ()
  return $ sqrt $ params V.! 1

constraint :: ScalarFunction (Double, Double)
constraint params g (a, b) = do
  case g of
    Just grad -> V.copy grad $ V.fromList [3 * a * t * t, -1.0]
    Nothing -> return ()
  return f
  where
    [x0, x1] = V.toList params
    t = (a * x0 + b)
    f = t * t * t - x1

lowerBounds :: V.Vector Double
lowerBounds = V.fromList [-1e222, 0.0]

checkReturn :: Result -> IO ()
checkReturn c = when (not $ isSuccess c) $ do
  putStrLn $ "NLOPT error '" ++ show c ++ "'!"
  E.exitFailure

xinit :: V.Vector Double
xinit = V.fromList [1.234, 5.678]

expectedParams :: V.Vector Double
expectedParams = V.fromList [0.33333337491594595,0.2962962707981303]

expectedMinimum :: Double
expectedMinimum = 0.5443310305302559

nnlsTest :: IO ()
nnlsTest = do
  putStrLn "NNLS"
  let a = V.fromList [1, 0, 1, 1, 1, 0, 0, 1, 1]
      b = V.fromList [1, 2, -1]
  (x, rNorm, mode) <- nnls 3 3 a b
  putStrLn $ "mode=" <> show mode
  putStrLn $ "x=" <> show x
  putStrLn $ "||Ax - b||_2 = " <> show rNorm


ldpTest :: IO ()
ldpTest = do
  putStrLn "LDP"
  let g = V.fromList [1, 0, 1, 1, 1, 0, 0, 1, 1]
      h = V.fromList [1, 1, 1]
  (x, xnorm, mode) <- ldp 3 3 g h
  putStrLn $ "mode=" <> show mode
  putStrLn $ "x=" <> show x
  putStrLn $ "||x||_2 = " <>  show xnorm

main :: IO ()
main = do
  nnlsTest
  ldpTest
  Just opt <- create LD_MMA 2
  putStrLn "Got optimizer."
  bnd <- set_lower_bounds opt lowerBounds
  checkReturn bnd
  putStrLn "Set lower bounds."
  obj <- set_min_objective opt objective ()
  checkReturn obj
  putStrLn "Set min objective."
  cnst1 <- add_inequality_constraint opt constraint (2, 0) 1e-8
  checkReturn cnst1
  putStrLn "Added inequality constraint #1."
  cnst2 <- add_inequality_constraint opt constraint (-1, 1) 1e-8
  checkReturn cnst2
  putStrLn "Added inequality constraint #2."
  tol <- set_xtol_rel opt 1e-4
  checkReturn tol
  putStrLn "Set param tolerance."
  Output result cost params <- optimize opt xinit
  checkReturn result
  let
    verr = V.sum . V.map (abs) $ V.zipWith (-) params expectedParams
    ferr = abs $ cost - expectedMinimum
  if verr > 1e-12 || ferr > 1e-12
    then do
    putStrLn "Solution error!"
    putStrLn $ "Expected minimum of " ++ show expectedMinimum ++
      " and got " ++ show cost
    putStrLn $ "Expected " ++ show expectedParams ++ " and got " ++ show params
    E.exitFailure
    else
    putStrLn $ "Found minimum at f" ++ show params ++ " = " ++ show cost
