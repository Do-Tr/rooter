#ifndef TipTimes_H
#define TipTimes_H

#include <map>
#include <string>
class Alignment;
class Settings;

struct Date {

    int d, m, y;
    Date(int a, int b, int c);
};



class TipTimes {

    public:
                                        TipTimes(void) = delete;
                                        TipTimes(Alignment* a, Settings* s);
        std::map<std::string, double>&  getTipTimesMap(void) { return tipTimesMap; }

    protected:
        int                             countLeapYears(Date d);
        int                             getDifference(Date dt1, Date dt2);
        bool                            isTaxonPresent(std::string tn);
        void                            parseTipDatesFile(std::string fn);
        Date                            processDateString(std::string t);
        std::map<std::string, double>   tipTimesMap;
};

#endif
