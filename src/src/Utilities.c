#include "config.h"
#include <sys/stat.h>

/*!
 * This function returns true (false) when a file exists (does not exist).
 */
bool CheckFile(const char FileName[]){
    FILE *fp;
    fp = fopen(FileName,"r");
    if(fp == NULL){
        return false;
    }
    fclose(fp);
    return true;
}

/*!
 * This function returns true (false) when a directory exists (does not exist).
 */
bool CheckDir(const char DirName[]){
    FILE *fp;
    fp = fopen(DirName,"r");
    if(fp == NULL){
        return false;
    }
    fclose(fp);
    return true;
}

/*!
 * This function generates a new directory. If the desired directory is found,
 * this function does not change the original directory. 
 */
void MakeDir(const char DirName[]){

    if(CheckDir(DirName) == false){
        int checkflag = mkdir(DirName,0755);
        if(checkflag < 0){
            fprintf(stderr,"Directory [%s] creation error.\n",DirName);
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr,"Directory [%s] is created.\n",DirName);
        }
    }

    return ;
}
