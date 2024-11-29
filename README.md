## What is alignshort?
alignshort is a very fast DNA short sequence alignment software that supports setting the number of mismatches and bidirectional alignment of nucleotides.

## Installation and usage
### Direct download
alignshort can run on windows and linux and is a static compiled program, can be downloaded directly: [Linux](https://raw.githubusercontent.com/gagayuan/alignshort/main/linux/alignshort), [Windows](https://raw.githubusercontent.com/gagayuan/alignshort/main/windows/alignshort.exe).
### Build from source
To build alignshort from source in Linux, you need gcc >= 11.4.0, make >= 4.3, and the commands for static compilation are as follows:
```
cd linux
make
make all
```
On Windows, the following is the process of statically compiling a 64-bit program based on the Visual Studio 2022 environment:
```
cd windows
& "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.39.33519\bin\Hostx64\x64\cl.exe" /MT ../alignshort.c ..\include\getopt9\src\getopt.c `
/I"..\include\pthreads_x64-windows-static\include" `
/I"..\include\getopt9\include" `
/I"..\include\msvc2017_64\include\zlib" `
/link /OUT:alignshort.exe `
/LIBPATH:"..\include\pthreads_x64-windows-static\lib" `
/LIBPATH:"..\include\msvc2017_64\lib\zlib" `
/LIBPATH:"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.39.33519\lib\x64" `
/LIBPATH:"C:\Program Files (x86)\Windows Kits\10\lib\10.0.22621.0\ucrt\x64" `
/LIBPATH:"C:\Program Files (x86)\Windows Kits\10\Lib\10.0.22621.0\um\x64" `
pthreadVC3.lib zlib.lib
```

### Using alignshort
The length of all short sequences for alignshort query fasta must be the same and less than 32. The query sequence must be in fasta format, and the target reference sequence file can be in fasta or fasta.gz format. The usage is as follows:
```
Usage: ./alignshort
        -q <query fasta file, Required parameters>
        -t <target fasta or fasta.gz file, Required parameters>
        -o <output file, in csv format, Required parameters>
        -m <mismatch number in alignment, Default:5>
        -n <number of threads, Default: all>
        -r <yes or no to use reverse complementary sequences of querySeq for alignment, Default: yes>

```

#### Result file format
The result file is in csv format and the columns are explained below:

| column name    | meaning |
|----------------|---------------------------------------|
| queryId        | ID of the query sequence |
| queryStrand    | direction of the sequence, `ori` means align with source sequence, `rc` means align with reverse complementary sequence |
| querySeq       | query sequence. If the queryIds are the same but the queryStrands are different, just means they are reverse complementary |
| matchId        | The ID of the target reference sequence to be aligned, only the characters before the first space of the reference ID are preserved.   |
| startPosition  | alignment start position on the reference sequence |
| endPosition    | alignment end position on the reference sequence |
| matchSeq       | matched sequence in reference sequence |
| mismatchNumber | number of nucleotide mismatches in the match |

### cpu and memory consumption

The performance of 1200 sequences (19nt), 5 mismatches bi-directional alignment to human cdna file(Homo_sapiens.GRCh38.cdna.all.fa, 429MB) was tested as follows:

| thread number | RunTime(second) | VSZ Memory(GB) | RSS Memory(GB) |
|---------------|-----------------|----------------|----------------|
| 255           | 45              | 91.9           | 14.1           |
| 200           | 48              | 74.4           | 14.1           |
| 125           | 51              | 51.0           | 14.1           |
| 100           | 51              | 43.1           | 14.1           |
| 64            | 72              | 31.9           | 14.1           |
| 32            | 114             | 23.8           | 14.1           |
| 16            | 222             | 18.8           | 14.1           |
| 8             | 465             | 16.3           | 14.1           |
| 4             | 1074            | 15.2           | 14.1           |

The performance of 1200 sequences (19nt), 5 mismatches bi-directional alignment to human dna file(Homo_sapiens.GRCh38.dna.toplevel.fa.gz, 3.2GB) was tested as follows:

| thread number | RunTime(second) | VSZ Memory(GB) | RSS Memory(GB) |
|---------------|-----------------|----------------|----------------|
| 255           | 618             | 176.2          | 108.4          |
| 200           | 603             | 168.2          | 108.4          |
| 125           | 591             | 142.6          | 108.4          |
| 100           | 600             | 136.9          | 108.4          |
| 64            | 753             | 125.9          | 108.4          |
| 32            | 1260            | 117.2          | 108.4          |
| 16            | 2256            | 112.9          | 108.4          |
| 8             | 4248            | 110.7          | 108.4          |
| 4             | 8952            | 109.6          | 108.4          |

Therefore, when align to human cdna, it is recommended to use a server with 16GB of memory.
When align to the human dna, it is recommended to use a server with 128GB of memory.

### Contact
If you have any questions, please contact <xuxi@iict.ac.cn>.