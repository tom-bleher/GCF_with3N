//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool 
*/
//____________________________________________________________________________

#include <iostream>
#include <vector>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string_view>
#include <unordered_map>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "log4cpp/SimpleLayout.hh"

#include <TSystem.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

using std::setw;
using std::setfill;
using std::cout;
using std::endl;
using std::vector;

using namespace genie;

bool genie::gAbortingInErr = false;

//____________________________________________________________________________
Messenger * Messenger::fInstance = 0;
//____________________________________________________________________________
Messenger::Messenger()
{
  fInstance =  0;
}
//____________________________________________________________________________
Messenger::~Messenger()
{
  fInstance = 0;
}
//____________________________________________________________________________
Messenger * Messenger::Instance()
{
  if(fInstance == 0) {

    // the first thing that get's printed in a GENIE session is the banner
    utils::print::PrintBanner();

    static Messenger::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new Messenger;

    log4cpp::Appender * appender;
    appender = new log4cpp::OstreamAppender("default", &cout);
    const char* layoutenv = gSystem->Getenv("GMSGLAYOUT");
    std::string layoutstr = (layoutenv) ? string(layoutenv) : "BASIC";
    if ( layoutstr == "SIMPLE" )
      appender->setLayout(new log4cpp::SimpleLayout());
    else
      appender->setLayout(new log4cpp::BasicLayout());


    log4cpp::Category & MSG = log4cpp::Category::getRoot();

    MSG.setAdditivity(false);
    MSG.addAppender(appender);

    fInstance->Configure(); // set user-defined priority levels
  }
  return fInstance;
}
//____________________________________________________________________________
namespace {
  // Messenger::operator() category cache (see below)
  bool gCategoryCacheOn = false;
  void DisableCategoryCache(void) { gCategoryCacheOn = false; }
}
log4cpp::Category & Messenger::operator () (const char * stream)
{
  // Every LOG() statement calls this, including the many whose priority is
  // then filtered out. log4cpp::Category::getInstance() locks a mutex, builds
  // a std::string and searches a std::map<std::string, Category*> each time,
  // which made it the largest single cost of INTRANUKE cascades. Categories
  // live until log4cpp's HierarchyMaintainer is destroyed at exit, so cache
  // their addresses (per thread, so no locking):
  //  - L1: direct-mapped on the stream pointer (LOG() streams are almost
  //    always string literals), verified by comparing the name, so a reused
  //    pointer can never return the wrong category;
  //  - L2: hash map keyed by the name.
  // An atexit hook that runs before the HierarchyMaintainer destructor
  // switches the cache off, after which this falls back to the plain lookup.
  static const bool cache_init = [] {
    log4cpp::Category::getRoot();       // construct the maintainer first, so
    std::atexit(DisableCategoryCache);  // this hook runs before it is destroyed
    gCategoryCacheOn = true;
    return true;
  }();
  (void) cache_init;
  if (!gCategoryCacheOn) return log4cpp::Category::getInstance(stream);

  struct Slot { const char * key; const char * name; log4cpp::Category * cat; };
  static thread_local Slot l1[256];        // trivially destructible
  Slot & slot = l1[(reinterpret_cast<std::uintptr_t>(stream) >> 3) & 255];
  if (slot.key == stream && std::strcmp(slot.name, stream) == 0) {
    return *slot.cat;
  }

  typedef std::unordered_map<std::string_view, log4cpp::Category *> Cache_t;
  static thread_local Cache_t * l2 = 0;    // never freed: no TLS destructor
  if (!l2) l2 = new Cache_t;
  log4cpp::Category * cat = 0;
  Cache_t::const_iterator it = l2->find(std::string_view(stream));
  if (it != l2->end()) {
    cat = it->second;
  } else {
    cat = &log4cpp::Category::getInstance(stream);
    l2->emplace(std::string_view(cat->getName()), cat);
  }
  slot.key  = stream;
  slot.name = cat->getName().c_str();
  slot.cat  = cat;
  return *cat;
}
//____________________________________________________________________________
void Messenger::SetPriorityLevel(
   const char * stream, log4cpp::Priority::Value priority)
{
  log4cpp::Category & MSG = log4cpp::Category::getInstance(stream);

  MSG.setPriority(priority);
}
//____________________________________________________________________________
void Messenger::Configure(void)
{
// The Configure() method will look for priority level xml config files, read
// the priority levels and set them.
// By default, first it will look for the $GENIE/config/Messenger.xml file.
// Look at this file for the  XML schema.

  string msg_config_file = utils::xml::GetXMLFilePath("Messenger.xml");

  // parse & set the default priority levels
  bool ok = this->SetPrioritiesFromXmlFile(msg_config_file);
  if(!ok) {
    SLOG("Messenger", pERROR)
      << "Priority levels from: " << msg_config_file << " were not set!";
  }

}
//____________________________________________________________________________
bool Messenger::SetPrioritiesFromXmlFile(string filenames)
{
// Reads input XML config file and sets the priority levels.
// The input can be a colection of XML files, delimited with a ":".
// The full path for each file should be given.
// In case of multiple files, all files are read.
// The later each file is, the higher priority it has (if the same mesg
// stream is listed twice with conflicting priority, the last one is used.).
//

  if(filenames.size()==0) return false;

  vector<string> filename_vec = utils::str::Split(filenames, ":");
  if(filename_vec.size() == 0) return false;

  vector<string>::const_iterator filename_iter;
  for(filename_iter  = filename_vec.begin();
      filename_iter != filename_vec.end(); ++filename_iter)
  {
     // look in all known XML locations, just like primary file
     string filename = utils::xml::GetXMLFilePath(*filename_iter);

    SLOG("Messenger", pNOTICE)
      << "Reading msg stream priorities from XML file: " << filename;

    xmlDocPtr xml_doc = xmlParseFile(filename.c_str());
    if(xml_doc==NULL) {
       SLOG("Messenger", pERROR)
           << "XML file could not be parsed! [file: " << filename << "]";
       return false;
    }

    xmlNodePtr xml_root = xmlDocGetRootElement(xml_doc);
    if(xml_root==NULL) {
       SLOG("Messenger", pERROR)
           << "XML doc. has null root element! [file: " << filename << "]";
       xmlFreeDoc(xml_doc);
       return false;
    }

    if( xmlStrcmp(xml_root->name, (const xmlChar *) "messenger_config") ) {
       SLOG("Messenger", pERROR)
         << "XML doc. has invalid root element! [file: " << filename << "]";
       xmlFreeNode(xml_root);
       xmlFreeDoc(xml_doc);
       return false;
    }

    xmlNodePtr xml_msgp = xml_root->xmlChildrenNode; // <priority>'s

    // loop over all xml tree nodes that are children of the <spline> node
    while (xml_msgp != NULL) {

      // enter everytime you find a <priority> tag
      if( (!xmlStrcmp(xml_msgp->name, (const xmlChar *) "priority")) ) {

         string msgstream = utils::str::TrimSpaces(
                  utils::xml::GetAttribute(xml_msgp, "msgstream"));
         string priority =
                utils::xml::TrimSpacesClean( xmlNodeListGetString(
                               xml_doc, xml_msgp->xmlChildrenNode, 1));
         log4cpp::Priority::Value pv = this->PriorityFromString(priority);
         this->SetPriorityLevel(msgstream.c_str(), pv);
         SLOG("Messenger", pINFO)
                  << "Set priority level: " << setfill('.')
                          << setw(24) << msgstream << " --> " << priority;
      }
      xml_msgp = xml_msgp->next;
    }//xml_msgp != NULL

    //xmlFree(xml_msgp);
    xmlFreeNode(xml_msgp);
    xmlFreeDoc(xml_doc);
  }//filename_iter

  return true;
}
//____________________________________________________________________________
log4cpp::Priority::Value Messenger::PriorityFromString(string p)
{
  if ( p.find("DEBUG")  != string::npos ) return log4cpp::Priority::DEBUG;
  if ( p.find("INFO")   != string::npos ) return log4cpp::Priority::INFO;
  if ( p.find("NOTICE") != string::npos ) return log4cpp::Priority::NOTICE;
  if ( p.find("WARN")   != string::npos ) return log4cpp::Priority::WARN;
  if ( p.find("ERROR")  != string::npos ) return log4cpp::Priority::ERROR;
  if ( p.find("CRIT")   != string::npos ) return log4cpp::Priority::CRIT;
  if ( p.find("ALERT")  != string::npos ) return log4cpp::Priority::ALERT;
  if ( p.find("FATAL")  != string::npos ) return log4cpp::Priority::FATAL;

  SLOG("Messenger", pWARN)
                    << "Unknown priority = " << p << " - Setting to INFO";
  return log4cpp::Priority::INFO;
}
//____________________________________________________________________________
