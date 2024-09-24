{-# LANGUAGE TemplateHaskell #-}


module EVM.CheatsTH where

import EVM.ABI

import Data.ByteString.Char8 (pack)
import Data.Vector qualified as V

import Language.Haskell.TH
import Language.Haskell.TH.Syntax


liftByteString :: String -> Q Exp
liftByteString txt = AppE (VarE 'pack) <$> lift txt

liftAbiType :: AbiType -> Q Exp
liftAbiType AbiBoolType = [| AbiBool |]
liftAbiType (AbiUIntType n) = [| AbiUInt $(lift n) |]
liftAbiType (AbiIntType n) = [| AbiInt $(lift n) |]
liftAbiType AbiAddressType = [| AbiAddress |]
liftAbiType (AbiBytesType n) = [| AbiBytes $(lift n) |]
liftAbiType AbiStringType = [| AbiString |]
liftAbiType AbiBytesDynamicType = [| AbiBytesDynamic |]
liftAbiType _ = error "unimplemented"

envReadSingleCheat :: String -> Q Exp
envReadSingleCheat sigString = [|
    \wrap convert ->
        action $sigL $
            \sig input -> case decodeBuf [AbiStringType] input of
                CAbi [AbiString variable] -> let
                    varStr = toString variable
                    cont value = continueOnce $ do
                        either' (convert value) frameRevert $ \v -> do
                            frameReturn $ wrap v
                    in query (PleaseReadEnv varStr cont)
                _ -> vmError (BadCheatCode sig)
    |]
    where
        sigL = liftByteString sigString

envReadMultipleCheat :: String -> AbiType -> Q Exp
envReadMultipleCheat sigString arrType = [|
    \convert ->
        action $sigL $
            \sig input -> case decodeBuf [AbiStringType, AbiStringType] input of
                CAbi [AbiString variable, AbiString delimiter] -> let
                    (varStr, delimStr) = (toString variable, toString delimiter)
                    cont value = continueOnce $ do
                        let (errors, values) = partitionEithers $ map convert $ splitOn delimStr value
                        case errors of
                            [] -> do
                                let result = AbiTuple $ V.fromList [AbiArrayDynamic $arrTypeL $ V.fromList $ map $wrapL values]
                                frameReturn result
                            (e:_) -> frameRevert e
                    in query (PleaseReadEnv varStr cont)
                _ -> vmError (BadCheatCode sig)
    |]
    where
        sigL = liftByteString sigString
        arrTypeL = liftData arrType
        wrapL = liftAbiType arrType