{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE TypeFamilies #-}

module EVM.Stack 
  ( Stack(.., Empty, (:<|))
  , empty
  , push
  , pop
  , peek
  , size
  , EVM.Stack.filter
  , EVM.Stack.fromList
  , EVM.Stack.splitAt
  , EVM.Stack.toList
  ) where

import Data.Monoid ()
import Data.List qualified as List
import Data.Semigroup ()
import GHC.Exts (IsList(..))
import Optics.Core (AsEmpty(..), nearly, IxValue, Ixed(..), atraversalVL, AffineTraversalVL', (<&>), Index)

-- | A stack data structure with O(1) operations
data Stack a = Stack
  { stackItems :: [a]  -- ^ The items in the stack (head is top)
  , stackSize  :: !Int -- ^ Cache of the stack size
  } deriving (Eq, Show)

-- | Empty stack pattern
pattern Empty :: Stack a
pattern Empty = Stack [] 0

-- | Pattern for matching on a non-empty stack
infixr 5 :<|
pattern (:<|) :: a -> Stack a -> Stack a
pattern x :<| xs <- (viewTop -> Just (x, xs))
  where
    x :<| xs = push x xs

-- | Instance for Monoid, allowing stacks to be combined with <>
instance Monoid (Stack a) where
  -- The empty stack
  mempty = empty
  
-- | Instance for Semigroup, required for Monoid
instance Semigroup (Stack a) where
  -- Concatenate two stacks, maintaining the correct size
  (Stack xs sz1) <> (Stack ys sz2) = Stack (xs <> ys) (sz1 + sz2)

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

-- | Instance for Foldable, allowing operations like foldr, foldl, toList, etc.
instance Foldable Stack where
  -- foldr applies a function to each element and an accumulator, from right to left
  foldr f z (Stack xs _) = foldr f z xs
  {-# INLINE foldr #-}
  
  -- foldl applies a function to each element and an accumulator, from left to right
  foldl f z (Stack xs _) = foldl f z xs
  {-# INLINE foldl #-}
  
  -- foldMap applies a function to each element and combines the results using a Monoid
  foldMap f (Stack xs _) = foldMap f xs
  {-# INLINE foldMap #-}
  
  -- null checks if the stack is empty
  null (Stack xs _) = null xs
  {-# INLINE null #-}
  
  -- length returns the number of elements in the stack
  length = size
  {-# INLINE length #-}
  
  -- elem checks if an element is in the stack
  elem x (Stack xs _) = elem x xs
  {-# INLINE elem #-}
  
  -- maximum finds the maximum element
  maximum (Stack xs _) = maximum xs
  {-# INLINE maximum #-}
  
  -- minimum finds the minimum element
  minimum (Stack xs _) = minimum xs
  {-# INLINE minimum #-}
  
  -- sum sums all elements
  sum (Stack xs _) = sum xs
  {-# INLINE sum #-}
  
  -- product multiplies all elements
  product (Stack xs _) = product xs
  {-# INLINE product #-}

viewTop :: Stack a -> Maybe (a, Stack a)
viewTop (Stack [] _) = Nothing
viewTop (Stack (x:xs) n) = Just (x, Stack xs (n-1))

-- | Create an empty stack
empty :: Stack a
empty = Stack [] 0

-- | Push an item onto the stack
push :: a -> Stack a -> Stack a
push x (Stack xs n) = Stack (x:xs) (n+1)

-- | Pop an item from the stack
pop :: Stack a -> Maybe (a, Stack a)
pop = viewTop

-- | Peek at the top item without popping
peek :: Stack a -> Maybe a
peek (Stack [] _) = Nothing
peek (Stack (x:_) _) = Just x

-- | Get the size of the stack
size :: Stack a -> Int
size s = s.stackSize

-- | Convert a list to a stack
fromList :: [a] -> Stack a
fromList xs = Stack xs (length xs)

-- | Convert a stack to a list
toList :: Stack a -> [a]
toList s = s.stackItems

-- | Split a stack at the specified index
-- Returns a tuple of two stacks: the first contains the first n elements,
-- the second contains the rest
splitAt :: Int -> Stack a -> (Stack a, Stack a)
splitAt n (Stack xs sz)
  | n <= 0    = (empty, Stack xs sz)
  | n >= sz   = (Stack xs sz, empty)
  | otherwise = (Stack front n, Stack back (sz - n))
  where
    (front, back) = List.splitAt n xs

-- | Filter elements in a stack based on a predicate
filter :: (a -> Bool) -> Stack a -> Stack a
filter p (Stack xs _) = 
  let filteredList = List.filter p xs
  in Stack filteredList (length filteredList)

-- Enable OverloadedLists extension support
instance IsList (Stack a) where
  type Item (Stack a) = a
  
  fromList = EVM.Stack.fromList
  toList = EVM.Stack.toList

-- | Instance for Functor, allowing mapping functions over stack elements
instance Functor Stack where
  -- fmap applies a function to each element in the stack
  fmap f (Stack xs n) = Stack (fmap f xs) n
  {-# INLINE fmap #-}
  
  -- (<$) replaces all elements with a constant value
  (<$) a (Stack xs n) = Stack (a <$ xs) n
  {-# INLINE (<$) #-}