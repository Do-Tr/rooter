#include <fstream>
#include <iostream>
#include <sstream>
#include "Alignment.h"
#include "Msg.h"
#include "Settings.h"
#include "TipTimes.h"

const int monthDays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};



TipTimes::TipTimes(Alignment* a, Settings* s) {

    // initialize the tip times from the alignment
    std::vector<std::string> taxa = a->getTaxonNames();
    for (std::string n : taxa)
        tipTimesMap.insert( std::make_pair(n, -1.0) );
    
    // parse the file containing the tip dates
    parseTipDatesFile(s->getTipTimesFileName());
}

bool TipTimes::isTaxonPresent(std::string tn) {

    std::map<std::string, double>::iterator it = tipTimesMap.find(tn);
    if (it == tipTimesMap.end())
        return false;
    return true;
}

void TipTimes::parseTipDatesFile(std::string fn) {

    // open the file
    std::ifstream ttStream(fn.c_str());
    if (ttStream.is_open() == true)
        std::cout << "   * Initializing the tip dates from file \"" << fn << "\"" << std::endl;
    else
        {
        std::cerr << "Cannot open file \"" + fn + "\"" << std::endl;
        exit(1);
        }

    std::string linestring = "";
    int line = 0;
    while ( getline (ttStream, linestring) )
        {
        std::istringstream linestream(linestring);
        //std::cout << line << " -- \"" << linestring << "\"" << std::endl;
        int ch;
        std::string word = "";
        int wordNum = 0;
        std::string tName = "";
        do
            {
            word = "";
            linestream >> word;
            wordNum++;
            //std::cout << "word: " << wordNum << "\" " << word << "\"" << std::endl;
            
            if (wordNum == 1)
                {
                // process taxon name
                if (isTaxonPresent(word) == false)
                    Msg::error("Taxon " + word + " in tip dates file not found");
                tName = word;
                }
            else if (wordNum == 2)
                {
                // process tip time
                Date b(18, 4, 1966);
                Date d = processDateString(word);
                int diff = getDifference(b, d);
                std::map<std::string, double>::iterator it = tipTimesMap.find(tName);
                if (it == tipTimesMap.end())
                    Msg::error("Taxon " + tName + " in tip dates map not found");
                it->second = (double)diff;
                }
            else
                {
                // shouldn't have three words on a line
                Msg::error("Found three columns on line" + std::to_string(line+1));
                }
            
            } while ( (ch=linestream.get()) != EOF );
            
        line++;
        }
    
    // close the file
    ttStream.close();
    
    // check that all of the taxa have a date
    for (std::map<std::string, double>::iterator it=tipTimesMap.begin(); it!=tipTimesMap.end(); it++)
        {
        if (it->second < 0.0)
            Msg::error("The sample date for " + it->first + " was never initialized");
        }
    
    // rescale all of the dates to the most recent, which has a time of 0.0
    double maxTime = 0.0;
    for (std::map<std::string, double>::iterator it=tipTimesMap.begin(); it!=tipTimesMap.end(); it++)
        {
        if (it->second > maxTime)
            maxTime = it->second;
        }
    for (std::map<std::string, double>::iterator it=tipTimesMap.begin(); it!=tipTimesMap.end(); it++)
        {
        it->second = maxTime - it->second;
        }
    
#   if 0
    for (std::map<std::string, double>::iterator it=tipTimesMap.begin(); it!=tipTimesMap.end(); it++)
        {
        std::cout << it->first << " " << it->second << std::endl;
        }
#   endif
}

int TipTimes::countLeapYears(Date d) {

    int years = d.y;
  
    // check if the current year needs to be considered for the count of leap years or not
    if (d.m <= 2)
        years--;
  
    // a year is a leap year if it is a multiple of 4, multiple of 400 and not a multiple of 100.
    return years / 4 - years / 100 + years / 400;
}

int TipTimes::getDifference(Date dt1, Date dt2) {
  
    // initialize count using years and day
    int n1 = dt1.y*365 + dt1.d;
  
    // Add days for months in given date
    for (int i=0; i<dt1.m - 1; i++)
        n1 += monthDays[i];
  
    // Since every leap year is of 366 days,
    // Add a day for every leap year
    n1 += countLeapYears(dt1);
  
    // SIMILARLY, COUNT TOTAL NUMBER OF DAYS BEFORE 'dt2'
  
    int n2 = dt2.y * 365 + dt2.d;
    for (int i=0; i<dt2.m - 1; i++)
        n2 += monthDays[i];
    n2 += countLeapYears(dt2);
  
    // return difference between two counts
    return (n2 - n1);
}

Date TipTimes::processDateString(std::string t) {

    Date date(0,0,0);
    
    // we assuem the string is year-month-day
    std::string s = "";
    int cnt = 0;
    for (int i=0; i<t.length(); i++)
        {
        char c = t[i];
        if (c == '-')
            {
            if (cnt == 0)
                date.y = std::stod(s);
            else if (cnt == 1)
                date.m = std::stod(s);
            s = "";
            cnt++;
            }
        if ( (s == "" && c == '0') || c == '-' )
            ;
        else
            s += c;
        }
    if (s != "")
        date.d = std::stod(s);

    //std::cout << date.m << " " << date.d << " " << date.y << std::endl;
    return date;
}

Date::Date(int a, int b, int c) {

    d = a;
    m = b;
    y = c;
}
