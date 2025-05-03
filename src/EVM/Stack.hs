{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ViewPatterns #-}

module EVM.Stack 
  ( Stack(.., Empty, (:<|))
  , push
  , size
  , EVM.Stack.filter
  , EVM.Stack.splitAt
  , EVM.Stack.toList
  ) where

import Data.Monoid ()
import Data.List qualified as List
import Data.Semigroup ()
import Optics.Core (AsEmpty(..), nearly, Index, IxValue, Ixed(..), atraversalVL, AffineTraversalVL', (<&>))

-- | A stack-like data structure with O(1) length, push, pop operations
data Stack a = Stack
  { stackItems :: [a]  -- ^ Items in the stack (head is top)
  , stackSize  :: !Int -- ^ Stack size
  } deriving (Eq, Show)

-- | Stack pattern synonyms

pattern Empty :: Stack a
pattern Empty = Stack [] 0

infixr 5 :<|
pattern (:<|) :: a -> Stack a -> Stack a
pattern x :<| xs <- (viewTop -> Just (x, xs))
  where
    x :<| xs = push x xs

viewTop :: Stack a -> Maybe (a, Stack a)
viewTop (Stack [] _) = Nothing
viewTop (Stack (x:xs) n) = Just (x, Stack xs (n-1))

-- | Stack operations

push :: a -> Stack a -> Stack a
push x (Stack xs n) = Stack (x:xs) (n+1)

size :: Stack a -> Int
size s = s.stackSize

toList :: Stack a -> [a]
toList s = s.stackItems

splitAt :: Int -> Stack a -> (Stack a, Stack a)
splitAt n (Stack xs sz)
  | n <= 0    = (Empty, Stack xs sz)
  | n >= sz   = (Stack xs sz, Empty)
  | otherwise = (Stack front n, Stack back (sz - n))
  where
    (front, back) = List.splitAt n xs

filter :: (a -> Bool) -> Stack a -> Stack a
filter p (Stack xs _) = 
  let filteredList = List.filter p xs
  in Stack filteredList (length filteredList)


-- | Base instances
instance Monoid (Stack a) where
  mempty = Empty

instance Semigroup (Stack a) where
  (Stack xs sz1) <> (Stack ys sz2) = Stack (xs <> ys) (sz1 + sz2)

instance Foldable Stack where
  elem x (Stack xs _) = List.elem x xs
  {-# INLINE elem #-}
  foldl f z (Stack xs _) = List.foldl f z xs
  {-# INLINE foldl #-}
  foldl' f z (Stack xs _) = List.foldl' f z xs
  {-# INLINE foldl' #-}
  foldl1 f (Stack xs _) = List.foldl1 f xs
  {-# INLINE foldl1 #-}
  foldr f z (Stack xs _) = List.foldr f z xs
  {-# INLINE foldr #-}
  foldr1 f (Stack xs _) = List.foldr1 f xs
  {-# INLINE foldr1 #-}
  foldMap f (Stack xs _) = foldMap f xs
  {-# INLINE foldMap #-}
  length = size
  {-# INLINE length #-}
  maximum (Stack xs _) = List.maximum xs
  {-# INLINE maximum #-}
  minimum (Stack xs _) = List.minimum xs
  {-# INLINE minimum #-}
  null (Stack _ n) = n == 0
  {-# INLINE null #-}

instance Functor Stack where
  fmap f (Stack xs n) = Stack (fmap f xs) n
  {-# INLINE fmap #-}
  a <$ (Stack xs n) = Stack (a <$ xs) n
  {-# INLINE (<$) #-}

-- | Optics instances
instance AsEmpty (Stack a) where
  _Empty = nearly Empty (\s -> size s == 0)
  {-# INLINE _Empty #-}

type instance Index (Stack a) = Int
type instance IxValue (Stack a) = a
instance Ixed (Stack a) where
  ix k = atraversalVL (ixListVL k)
  {-# INLINE ix #-}

ixListVL :: Int -> AffineTraversalVL' (Stack a) a
ixListVL k point f s@(Stack xs0 l) =
  if k < 0
  then point s
  else let go [] _ = point []
           go (a:as) 0 = f a <&> (:as)
           go (a:as) i = (a:) <$> (go as $! i - 1)
       in  (`Stack` l) <$> (go xs0 k)
{-# INLINE ixListVL #-}