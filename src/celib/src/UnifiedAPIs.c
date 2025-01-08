#include "config.h"

/*!
 * This function retunrs the feedback results. 
 */
struct CELibStructFeedbackOutput CELibGetFeedback(struct CELibStructFeedbackInput Input, const int Type){

    switch (Type){
        case CELibFeedbackType_SNII:
            return CELibGetSNIIFeedback((struct CELibStructFeedbackInput){
                    .Mass = Input.Mass,
                    .Metallicity = Input.Metallicity,
                    .MassConversionFactor = Input.MassConversionFactor,
                    .Elements = Input.Elements,
                    .noPopIII = Input.noPopIII,
                    });
            break;
        case CELibFeedbackType_SNIa:
            return CELibGetSNIaFeedback((struct CELibStructFeedbackInput){
                    .Mass = Input.Mass,
                    .Metallicity = Input.Metallicity,
                    .MassConversionFactor = Input.MassConversionFactor,
                    });
            break;
        case CELibFeedbackType_AGB:
            return CELibGetAGBFeedback((struct CELibStructFeedbackInput){
                    .Mass = Input.Mass,
                    .Metallicity = Input.Metallicity,
                    .MassConversionFactor = Input.MassConversionFactor,
                    .Elements = Input.Elements,
                    .Count = Input.Count,
                    .noPopIII = Input.noPopIII,
                    });
            break;
        case CELibFeedbackType_NSM:
            return CELibGetNSMFeedback((struct CELibStructFeedbackInput){
                    .Mass = Input.Mass,
                    .MassConversionFactor = Input.MassConversionFactor,
                    });
            break;
        default:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            break;
    }

    return (struct CELibStructFeedbackOutput){
        .Energy = 0,
        .EjectaMass = 0,
        .RemnantMass = Input.Mass*Input.MassConversionFactor,
        .Elements[0] = 0,
        .Elements[1] = 0,
        .Elements[2] = 0,
        .Elements[3] = 0,
        .Elements[4] = 0,
        .Elements[5] = 0,
        .Elements[6] = 0,
        .Elements[7] = 0,
        .Elements[8] = 0,
        .Elements[9] = 0,
        .Elements[10] = 0,
        .Elements[11] = 0,
        .Elements[12] = 0,
        };
}


/*!
 * This function returns next event time. 
 */
double CELibGetNextEventTime(struct CELibStructNextEventTimeInput Input, const int Type){

    switch (Type){
        case CELibFeedbackType_SNII:
            {
                double Time = CELibGetSNIIExplosionTime(Input.R,Input.Metallicity);
                if(Time < 0.0){
                    fprintf(stderr,"Input rate is incorrect.\n");
                    return 0.e0;
                } else { 
                    return Time;
                }
            }
        case CELibFeedbackType_SNIa:
            return CELibGetSNIaExplosionTime(Input.R,Input.Metallicity,Input.InitialMass_in_Msun,Input.Count);   
        case CELibFeedbackType_AGB:
            return CELibGetAGBFeedbackTime(Input.R,Input.Metallicity,Input.Count);   
        case CELibFeedbackType_NSM:
            return CELibGetNSMFeedbackTime(Input.R,Input.Metallicity,Input.InitialMass_in_Msun,Input.Count);   
        default:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            return 0.e0;
    }
}


/*!
 * This function retunrs the star-by-star feedback results. 
 */
struct CELibStructFeedbackStarbyStarOutput CELibGetFeedbackStarbyStar(struct CELibStructFeedbackStarbyStarInput Input, const int Type){

    switch (Type){
        case CELibFeedbackType_SNII:
            return CELibGetSNIIFeedbackStarbyStar((struct CELibStructFeedbackStarbyStarInput){
                    .Mass = Input.Mass,
                    .Metallicity = Input.Metallicity,
                    .MassConversionFactor = Input.MassConversionFactor,
                    .Elements = Input.Elements,
                    .noPopIII = Input.noPopIII,
                    });
            break;
        case CELibFeedbackType_SNIa:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            break;
        case CELibFeedbackType_AGB:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            break;
        case CELibFeedbackType_NSM:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            break;
        default:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            break;
    }

    //return ;
    return (struct CELibStructFeedbackStarbyStarOutput){
        .Energy = 0,
        .EjectaMass = 0,
        .RemnantMass = Input.Mass*Input.MassConversionFactor,
        .Elements[0] = 0,
        .Elements[1] = 0,
        .Elements[2] = 0,
        .Elements[3] = 0,
        .Elements[4] = 0,
        .Elements[5] = 0,
        .Elements[6] = 0,
        .Elements[7] = 0,
        .Elements[8] = 0,
        .Elements[9] = 0,
        .Elements[10] = 0,
        .Elements[11] = 0,
        .Elements[12] = 0,
        };
}


/*!
 * This function returns next event time. 
 */
double CELibGetNextEventTimeStarbyStar(struct CELibStructNextEventTimeStarbyStarInput Input, const int Type){

    switch (Type){
        case CELibFeedbackType_SNII:
            {
                /*
                double Time;
                if(Input.noPopIII == 1){
                    Time = CELibGetLifeTimeofStar(Input.InitialMass_in_Msun,CELibLifeTimeZ[1]);
                } else {
                    Time = CELibGetLifeTimeofStar(Input.InitialMass_in_Msun,Input.Metallicity);
                }
                */
                double Time = CELibGetLifeTimeofStar(Input.InitialMass_in_Msun,Input.Metallicity);
                if(Time < 0.0){
                    fprintf(stderr,"Input rate is incorrect.\n");
                    return 0.e0;
                } else { 
                    return Time;
                }
            }
        case CELibFeedbackType_SNIa:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            return 0.e0;
        case CELibFeedbackType_AGB:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            return 0.e0;
        case CELibFeedbackType_NSM:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            return 0.e0;
        default:
            fprintf(stderr,"Incorrect feedback type is used.\n");
            return 0.e0;
    }
}
