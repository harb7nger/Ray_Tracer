The program generates a pattern of RGB colors by taking weighted values of RGB.

How to run the program:

- the program can be compiled as: 
$ g++ hw0.cpp

- the input file is of the format:
$ cat input.txt
imsize 640 480

- the output ppm file can be generated as:
$ ./a.exe input.txt

- the output file would be in ppm format, in the same folder and can be opened using gimp
"inputfilename".ppm

What the program does:

- checks if the input is valid or not

- gets the dimensions of the image from the input file
    - parses the input file to get the dimensions 
    - returns invalid values if dimensions or input are not as expected
         
 
- creates output file and writes the header information

- drawsPatterns to the output file of the specified dimensions
    - calls a helper function to generate RGB pattern, the RGB pattern has weighted values of RGB
    - writes it to the output file which is "inputFileName".ppm
 
    