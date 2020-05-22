#ifndef POINTSOURCE_INCLUDE_TRACKSTAT_H_
#define POINTSOURCE_INCLUDE_TRACKSTAT_H_
#include "TString.h"
#include <string.h>

extern double dRA, dDEC;

class ArgInit {
public:
  TString trackstatroot;
  int source_type;
  ArgInit(int argc, char **argv) {
    if (argc < 3) {
      printf("%s  trackstat.root  SourceType  [source RA]  [source "
             "DEC]\n(SourceType: -m for Moon, -s for Sun, -p need RA & DEC)\n",
             argv[0]);
      exit(0);
    }
    trackstatroot = argv[1];
    source_type = 0;
    for (int i = 2; i < argc; ++i) {
      if (strcmp(argv[i], "-m") == 0) {
        source_type = 1;
      } else if (strcmp(argv[i], "-s") == 0) {
        source_type = 2;
      } else if (strcmp(argv[i], "-p") == 0) {
        source_type = 3;
        if (i + 3 <= argc) {
          dRA = atof(argv[++i]);
          dDEC = atof(argv[++i]);
        } else {
          printf("wrong RA, DEC!\n");
          exit(0);
        }
      }
    }
    if (!source_type) {
      printf("source error!\n");
      exit(0);
    }
  };
};
#endif
