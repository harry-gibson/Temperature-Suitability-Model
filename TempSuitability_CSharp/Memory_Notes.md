Approximately 700 files
Files are 32-bit float but we will probably read them into double arrays
5600 bytes per pixel

Work with 256 pixel square tiles
= 350 Mb per tile

We don't need to store the daily or 2 hourly temperature estimates in memory: - 
we only need each once, so there is no advantage in doing so. 

Therefore we run the model for the entire period in one pass! 
