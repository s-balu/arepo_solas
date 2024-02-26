struct CELibStructRunParameters CELibRunParameters = {
    .IntegrationSteps = 10000,
    .LifeTimeType = 1,

    /* Parameters for SNeII */
    .SNIIYieldsTableID = CELibSNIIYieldsTableID_P98
};

void CELibInit(void){

    CELibInitLifeTime();

    CELibInitSNIIYields();

    return ;
}