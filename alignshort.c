#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <pthread.h>
#include <getopt.h>
#include <stdbool.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#pragma warning(disable:4996)


#ifdef _WIN32
#define strdup _strdup
#include <windows.h>
int getNumberOfCores() {
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
}
#elif __linux__
#include <unistd.h>
int getNumberOfCores() {
    return sysconf(_SC_NPROCESSORS_ONLN);
}
#else
#error "Unsupported OS"
#endif

typedef struct {
    char* name;
    uint64_t sirnaSeqCodedInt;
} SIRNA;

typedef struct {
    char* name;
    char* sequence;
    uint32_t indexNum;
} Gene;

typedef struct {
    char* name;
} Gene2;

typedef struct {
    uint64_t shortSeqCodedInt;
    uint32_t idNumber;
    uint32_t idStart;
} ShortSeqStruct;

char* reverseComplement(const char* dna) {
    int length = strlen(dna);
    char* reversed = malloc(length + 1);
    if (reversed == NULL) {
        return NULL;
    }
    reversed[length] = '\0';
    for (int i = 0; i < length; i++) {
        switch (dna[length - 1 - i]) {
        case 'A':
            reversed[i] = 'T';
            break;
        case 'T':
            reversed[i] = 'A';
            break;
        case 'C':
            reversed[i] = 'G';
            break;
        case 'G':
            reversed[i] = 'C';
            break;
        default:
            reversed[i] = 'N';
            break;
        }
    }
    return reversed;
}

// parse fasta or fasta.gz file
Gene* parseFasta(const char* filename, uint32_t* count) {
    gzFile fp;
    kseq_t* seq;
    Gene* genes = NULL;

    fp = gzopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Could not open fasta.gz file %s\n", filename);
        return NULL;
    }

    seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        Gene* temp = realloc(genes, (*count + 1) * sizeof(Gene));
        if (!temp) {
            fprintf(stderr, "Failed to allocate memory\n");
            free(genes);
            kseq_destroy(seq);
            gzclose(fp);
            return NULL;
        }
        genes = temp;

        genes[*count].name = strdup(seq->name.s);
        genes[*count].sequence = strdup(seq->seq.s);
        genes[*count].indexNum = *count;

        (*count)++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    return genes;
}

int checkEqualLengthFasta(const char* filename) {
    gzFile fp;
    kseq_t* seq;
    int firstLength = 0;
    int isEqual = 1;

    fp = gzopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Could not open fasta.gz file %s\n", filename);
        return 0;
    }

    seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        int currentLength = seq->seq.l;

        if (firstLength == 0) {
            firstLength = currentLength;
        }
        else if (currentLength != firstLength) {
            isEqual = 0;
            break;
        }
    }

    kseq_destroy(seq);
    gzclose(fp);

    return isEqual ? firstLength : 0;
}


void printCurrentTime() {
    time_t now;
    struct tm* tm_now;
    char buf[80];
    time(&now);
    tm_now = localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", tm_now);
    printf("[%s] ", buf);
}

//Sequence encoding function, the four characters of ACGT are encoded in binary 00, 01, 10, 11, respectively, 
//then the 19bp short sequence of cDNA can be converted to 38bit, which can be stored in a 64bit int
uint64_t encodeSequence(const char* sequence) {
    uint64_t coded = 0;
    int stepLength = strlen(sequence);

    for (int j = 0; j < stepLength; j++) {
        coded <<= 2;
        switch (sequence[j]) {
        case 'A': break;
        case 'C': coded |= 1; break;
        case 'G': coded |= 2; break;
        case 'T': coded |= 3; break;
        default:
            fprintf(stderr, "Invalid character %c in DNA sequence\n", sequence[j]);
            return 0;
        }
    }
    return coded;
}

char* decodeSequence(uint64_t coded, int seqLen, char* sequence) {
    sequence[seqLen] = '\0';
    for (int i = 0; i < seqLen; i++) {
        uint8_t bits = coded & 0x3;
        switch (bits) {
        case 0: sequence[seqLen - 1 - i] = 'A'; break;
        case 1: sequence[seqLen - 1 - i] = 'C'; break;
        case 2: sequence[seqLen - 1 - i] = 'G'; break;
        case 3: sequence[seqLen - 1 - i] = 'T'; break;
        }
        coded >>= 2;
    }
    return sequence;
}


SIRNA* convertGenesToSiRNA(Gene* genes, int numGenes, const char* rcmap) {
    if (strcmp(rcmap, "yes") == 0) {
        SIRNA* sirnas = malloc(2 * numGenes * sizeof(SIRNA));
        for (int i = 0; i < numGenes; i++) {
            int neededSize = strlen(genes[i].name) + strlen(",ori") + 1;
            sirnas[i].name = malloc(neededSize);
            if (sirnas[i].name != NULL) {
                strcpy(sirnas[i].name, genes[i].name);
                strcat(sirnas[i].name, ",ori");
            }
            sirnas[i].sirnaSeqCodedInt = encodeSequence(genes[i].sequence);
        }
        for (int i = numGenes; i < 2 * numGenes; i++) {
            int neededSize = strlen(genes[i - numGenes].name) + strlen(",rc") + 1;
            sirnas[i].name = malloc(neededSize);
            if (sirnas[i].name != NULL) {
                strcpy(sirnas[i].name, genes[i - numGenes].name);
                strcat(sirnas[i].name, ",rc");
            }
            sirnas[i].sirnaSeqCodedInt = encodeSequence(reverseComplement(genes[i - numGenes].sequence));
        }
        return sirnas;
    }
    else {
        SIRNA* sirnas = malloc(numGenes * sizeof(SIRNA));
        for (int i = 0; i < numGenes; i++) {
            int neededSize = strlen(genes[i].name) + strlen(",ori") + 1;
            sirnas[i].name = malloc(neededSize);
            if (sirnas[i].name != NULL) {
                strcpy(sirnas[i].name, genes[i].name);
                strcat(sirnas[i].name, ",ori");
            }
            sirnas[i].sirnaSeqCodedInt = encodeSequence(genes[i].sequence);
        }
        return sirnas;
    }
}

const uint64_t m1 = 0x5555555555555555;
const uint64_t m2 = 0x3333333333333333;
const uint64_t m4 = 0x0f0f0f0f0f0f0f0f;
const uint64_t m8 = 0x00ff00ff00ff00ff;
const uint64_t m16 = 0x0000ffff0000ffff;
const uint64_t m32 = 0x00000000ffffffff;
const uint64_t evenMask = 0xAAAAAAAAAAAAAAAAULL;
uint8_t popcount64a(uint64_t x)
{
    x = (x & m2) + ((x >> 2) & m2);
    x = (x & m4) + ((x >> 4) & m4);
    x = (x & m8) + ((x >> 8) & m8);
    x = (x & m16) + ((x >> 16) & m16);
    x = (x & m32) + ((x >> 32) & m32);
    return x;
}

//Bit-or operations on neighboring bits in 64-bit involve Hamming weight calculations 
uint8_t computeMismatch(uint64_t xorResult) {
    return popcount64a((xorResult & m1) | ((xorResult & evenMask) >> 1));
}

typedef struct {
    Gene2* genes;
    SIRNA* sirnaArray;
    int sirnaCount;
    int mismatchNum;
    int shortSeqLen;
    ShortSeqStruct* sequences;
    uint64_t start;
    uint64_t end;
    char*** results;
    uint64_t* resultCounts;
} ThreadData;

//thread Function for sequence alignment
void* processPart(void* arg) {
    ThreadData* data = (ThreadData*)arg;

    char* sequence = malloc(data->shortSeqLen + 1);
    char* query_sequence = malloc(data->shortSeqLen + 1);
    char* str_line = (char*)malloc(150);

    for (uint64_t i = data->start; i < data->end; i++) {
        ShortSeqStruct seq = data->sequences[i];

        data->results[i - data->start] = NULL;
        uint64_t results_count = 0;

        for (int j = 0; j < data->sirnaCount; j++) {
            uint64_t xorResult = seq.shortSeqCodedInt ^ data->sirnaArray[j].sirnaSeqCodedInt;
            uint8_t mismatches = computeMismatch(xorResult);

            if (mismatches <= data->mismatchNum) {
                decodeSequence(seq.shortSeqCodedInt, data->shortSeqLen, sequence);
                decodeSequence(data->sirnaArray[j].sirnaSeqCodedInt, data->shortSeqLen, query_sequence);
                sprintf(str_line, "%s,%s,%s,%d,%d,%s,%u", data->sirnaArray[j].name, query_sequence, data->genes[seq.idNumber].name, seq.idStart + 1, seq.idStart + data->shortSeqLen, sequence, mismatches);
                char** new_results = realloc(data->results[i - data->start], (results_count + 1) * sizeof(char*));
                if (new_results) {
                    data->results[i - data->start] = new_results;
                    data->results[i - data->start][results_count++] = strdup(str_line);
                }
            }
        }
        data->resultCounts[i - data->start] = results_count;
    }
    return NULL;
}

//thread create and join Function for sequence Alignment
void processShortSeqs(Gene2* genes, uint64_t numStructs, ShortSeqStruct* sequences, SIRNA* sirnaArray, int sirnaCount, int mismatchNum, const char* outputFile, int shortSeqLen, int numThreads) {
    uint64_t structsPerThread = numStructs / numThreads;
    int remainder = numStructs % numThreads;
    int actualNumThreads = numThreads + (remainder > 0 ? 1 : 0);
    pthread_t* threads = malloc(actualNumThreads * sizeof(pthread_t));
    ThreadData* threadData = malloc(actualNumThreads * sizeof(ThreadData));

    char*** allResults = malloc(numStructs * sizeof(char**));
    uint64_t* allResultsCount = calloc(numStructs, sizeof(uint64_t));

    for (uint64_t i = 0; i < numStructs; i++) {
        allResults[i] = NULL;
    }

    for (int i = 0; i < actualNumThreads; i++) {
        threadData[i].genes = genes;
        threadData[i].sirnaArray = sirnaArray;
        threadData[i].sirnaCount = sirnaCount;
        threadData[i].mismatchNum = mismatchNum;
        threadData[i].shortSeqLen = shortSeqLen;
        threadData[i].sequences = sequences;
        threadData[i].start = i * structsPerThread;
        threadData[i].end = (i + 1) * structsPerThread;
        if (i == actualNumThreads - 1 && remainder != 0) {
            threadData[i].end = numStructs;
        }
        threadData[i].results = allResults + threadData[i].start;
        threadData[i].resultCounts = allResultsCount + threadData[i].start;
        int create_status = pthread_create(&threads[i], NULL, processPart, &threadData[i]);
        if (create_status != 0) {
            perror("Thread creation failed");
            return;
        }
    }

    for (int i = 0; i < actualNumThreads; i++) {
        int join_status = pthread_join(threads[i], NULL);
        if (join_status != 0) {
            perror("Thread join failed");
        }
    }

    printCurrentTime();
    printf("Start write to out file \n");

    FILE* outFile = fopen(outputFile, "wb");
    char* headstr = "queryId,queryStrand,querySeq,matchId,startPosition,endPosition,matchSeq,mismatchNumber";
    fprintf(outFile, "%s\n", headstr);
    for (uint64_t i = 0; i < numStructs; i++) {
        if (allResults[i]) {
            if (allResultsCount[i] > 0) {
                for (uint64_t j = 0; j < allResultsCount[i]; j++) {
                    fprintf(outFile, "%s\n", allResults[i][j]);
                    free(allResults[i][j]);
                }
                free(allResults[i]);
            }
        }
    }
    fclose(outFile);
    free(sequences);
    free(allResults);
    free(allResultsCount);
    free(threads);
    free(threadData);
}

typedef struct {
    Gene* genes;
    uint32_t geneStart;
    uint32_t geneEnd;
    int stepLength;
    ShortSeqStruct** shortSeqs;
    int* numShortSeqs;
} ThreadData2;

//thread Function for Constructing(step and encode) short sequence databases
void seqWalkAndCode(uint32_t idNum, Gene gene, int stepLength, ShortSeqStruct** shortSeqs, int* numShortSeqs) {
    int sequenceLength = strlen(gene.sequence);
    int initial_numShortSeqs = sequenceLength - stepLength + 1;
    if (initial_numShortSeqs <= 0) return;
    *shortSeqs = realloc(*shortSeqs, sizeof(ShortSeqStruct) * initial_numShortSeqs);

    for (uint32_t i = 0; i < initial_numShortSeqs; i++) {
        uint64_t coded = 0;
        bool hasN = false;
        for (int j = 0; j < stepLength; j++) {
            coded <<= 2;
            switch (gene.sequence[i + j]) {
            case 'N': hasN = true; break;
            case 'A': break;
            case 'C': coded |= 1; break;
            case 'G': coded |= 2; break;
            case 'T': coded |= 3; break;
            }
        }
        if (hasN) continue;
        (*shortSeqs)[*numShortSeqs].idNumber = idNum;
        (*shortSeqs)[*numShortSeqs].idStart = i;
        (*shortSeqs)[*numShortSeqs].shortSeqCodedInt = coded;
        (*numShortSeqs)++;
    }
}

void* threadFunction(void* arg) {
    ThreadData2* data = (ThreadData2*)arg;
    for (uint32_t i = data->geneStart; i < data->geneEnd; i++) {
        data->shortSeqs[i] = NULL;
        data->numShortSeqs[i] = 0;
        seqWalkAndCode(i, data->genes[i], data->stepLength, &data->shortSeqs[i], &data->numShortSeqs[i]);
    }
    return NULL;
}

int main(int argc, char* argv[]) {
    int opt;
    char* fastaPath = NULL;
    char* queryfastaFilePath = NULL;
    char* outputFile = NULL;
    int mismatchNum = 5;
    int numThreads = getNumberOfCores() - 1;
    char* rcmap = "yes";

    const char* usage = "Usage: %s \n\
        -q <query fasta file, Required parameters>\n\
        -t <target fasta file, Required parameters>\n\
        -o <output file, in csv format, Required parameters>\n\
        -m <mismatch number in alignment, Default:5>\n\
        -n <number of threads, Default: all>\n\
        -r <yes or no to use reverse complementary sequences of querySeq for alignment, Default: yes>\n";

    while ((opt = getopt(argc, argv, "q:t:o:m:n:r:")) != -1) {
        switch (opt) {
        case 'q':
            queryfastaFilePath = optarg;
            break;
        case 't':
            fastaPath = optarg;
            break;
        case 'o':
            outputFile = optarg;
            break;
        case 'm':
            mismatchNum = atoi(optarg);
            break;
        case 'n':
            numThreads = atoi(optarg) - 1;
            break;
        case 'r':
            rcmap = optarg;
            break;
        default:
            fprintf(stderr, usage, argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (!queryfastaFilePath || !fastaPath || !outputFile) {
        fprintf(stderr, "Missing mandatory arguments.\n");
        fprintf(stderr, usage, argv[0]);
        return EXIT_FAILURE;
    }
    printf("Input arguments: \n");
    printf("Query File Path: %s\n", queryfastaFilePath);
    printf("Target File Path: %s\n", fastaPath);
    printf("Output File: %s\n", outputFile);
    printf("Mismatch Number: %d\n", mismatchNum);
    printf("Number of Threads: %d\n", numThreads);
    printf("Whether to use reverse complementary sequences of querySeq for alignment: %s\n\n", rcmap);

    //check all query seq whether have same length
    int stepLength = checkEqualLengthFasta(queryfastaFilePath);
    if (stepLength == 0) {
        printf("ERROR! The query fasta sequences have different lengths !\n");
        return EXIT_FAILURE;
    }
    else {
        printf("Detected all query fasta sequences have equal length: %d\n", stepLength);
    }


    printCurrentTime();
    printf("Start Align short seq \n");

    // Construct(step and encode) short sequence database
    uint32_t geneCount = 0;
    Gene* genes = parseFasta(fastaPath, &geneCount);
    if (genes == NULL) {
        fprintf(stderr, "Failed to parse fasta file.\n");
        return EXIT_FAILURE;
    }

    printCurrentTime();
    printf("Parsed target fasta %s\n", fastaPath);

    int* totalNumShortSeqs = malloc(geneCount * sizeof(int));
    ShortSeqStruct** globalShortSeqs = malloc(geneCount * sizeof(ShortSeqStruct*));

    pthread_t* threads = malloc(numThreads * sizeof(pthread_t));
    int genesPerThread = geneCount / numThreads;
    int remainingGenes = geneCount % numThreads;

    for (int i = 0; i < numThreads; i++) {
        ThreadData2* data = malloc(sizeof(ThreadData2));
        data->genes = genes;
        data->geneStart = i * genesPerThread;
        data->geneEnd = (i + 1) * genesPerThread + (i == numThreads - 1 ? remainingGenes : 0);
        data->stepLength = stepLength;
        data->shortSeqs = globalShortSeqs;
        data->numShortSeqs = totalNumShortSeqs;

        if (pthread_create(&threads[i], NULL, threadFunction, data) != 0) {
            perror("Failed to create thread");
            return 1;
        }
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);

    // flatten array
    uint64_t totalNumShortSeqs2 = 0;
    for (int i = 0; i < geneCount; i++) {
        totalNumShortSeqs2 += totalNumShortSeqs[i];
    }
    printCurrentTime();
    printf("totalNumShortSeqs: %lld\n", totalNumShortSeqs2);
    printCurrentTime();
    printf("totalNumShortSeqs needed mem: %lld bytes (%.3f MB)\n", totalNumShortSeqs2 * sizeof(ShortSeqStruct), (double)(totalNumShortSeqs2 * sizeof(ShortSeqStruct)) / 1048576);
    ShortSeqStruct* globalShortSeqs2 = malloc(totalNumShortSeqs2 * sizeof(ShortSeqStruct));
    if (!globalShortSeqs2) {
        perror("Failed to allocate memory for globalShortSeqs2");
        exit(EXIT_FAILURE);
    }

    uint64_t currentIndex = 0;
    for (int i = 0; i < geneCount; i++) {
        if (globalShortSeqs[i] != NULL) {
            memcpy(&globalShortSeqs2[currentIndex], globalShortSeqs[i], totalNumShortSeqs[i] * sizeof(ShortSeqStruct));
            currentIndex += totalNumShortSeqs[i];
        }
    }

    Gene2* genes2 = malloc(geneCount * sizeof(Gene2));
    for (int i = 0; i < geneCount; i++) {
        genes2[i].name = strdup(genes[i].name);
    }
    for (int i = 0; i < geneCount; i++) {
        free(genes[i].name);
        free(genes[i].sequence);
    }
    free(genes);
    for (int i = 0; i < geneCount; i++) {
        if (globalShortSeqs[i] != NULL) {
            free(globalShortSeqs[i]);
        }
    }

    free(totalNumShortSeqs);
    free(globalShortSeqs);

    printCurrentTime();
    printf("Generated all short sequence database with a length of %d\n", stepLength);

    // Align all query seq to short sequence database
    int sirnaCount = 0;
    Gene* ori_sirnas = parseFasta(queryfastaFilePath, &sirnaCount);
    if (ori_sirnas == NULL) {
        fprintf(stderr, "Failed to parse query fasta file.\n");
        return EXIT_FAILURE;
    }
    SIRNA* sirnas = convertGenesToSiRNA(ori_sirnas, sirnaCount, rcmap);
    if (strcmp(rcmap, "yes") == 0) sirnaCount *= 2;
    processShortSeqs(genes2, totalNumShortSeqs2, globalShortSeqs2, sirnas, sirnaCount, mismatchNum, outputFile, stepLength, numThreads);

    printCurrentTime();
    printf("Finish Align short seq \n");
    return EXIT_SUCCESS;
}