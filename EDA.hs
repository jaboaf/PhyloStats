import Data.List
import Data.List
import System.IO

-- Open the file using ReadMode
-- Get the contents of the file
-- Close the file
readfromfile = do
    theFile <- openFile "Data/allnuc0712.tar/allnuc0712/EPI_ISL_479737.fasta" ReadMode
    contents <- hGetContents theFile
    putStr contents
    hClose theFile
