// ======================================================================
//
// Name: DummySource_source.cc
//
// Purpose: Source module that returns one or more empty events per input file.
//          Input files are not actually opened.
//
// FCL parameters:
//
// MaxEvents      - Number of empty events to generate per input file (default 1).
// CopyInput      - Boolean (default false).  If true, copy input file to current directory.
// CopyErrorFatal - Boolean (default true) - If true treat copy input errors as fatal
//                  (see usage notes).
//
// Created: 27-Feb-2025  H. Greenlee
//
// Usage notes.
//
// 1.  The intended purpose of this source module is to make input file(s)
//     accessible to later processing, but not to do anything else with them.
//     It does this in two ways.
//
//     a) By storing the name of the input file(s) in the internal sam parentage metadata
//        of generated root output file (if input was read from sam and RootOutput is
//        configured).
//
//     b) By copying the input file to the current directory (if fcl parameter
//        CopyInput is true).
//
// 2.  Copy input errors are handled in the following way.
//
//     a) If input file is already in the current directory, don't try to
//        copy (not an error).
//
//     b) It is an error if the input file appears to be a url (e.g. xrootd url).
//
//     c) Any kind of error returned returned by the copy command is treated as an error.
//
//     In cases (b) and (c), an error message will be printed, and a c++ exception will
//     be thrown if fcl parameter CopyErrorFatal is true.
//
// ======================================================================

#include <filesystem>
#include <sstream>
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"

namespace ubcore {
  class DummySourceDetail
  {
  public:

    // Constructor.
    DummySourceDetail(fhicl::ParameterSet const&,
                     art::ProductRegistryHelper&,
                     art::SourceHelper const&);

    // Required methods.
    void readFile(std::string const& filename, art::FileBlock*& fb);
    bool readNext(art::RunPrincipal const* const inR,
                  art::SubRunPrincipal const* const inSR,
                  art::RunPrincipal*& outR,
                  art::SubRunPrincipal*& outSR,
                  art::EventPrincipal*& outE);
    void closeCurrentFile();

  private:

    // Fcl parameters

    unsigned int fMaxEvents;   // Number of empty events to return per input file.
    bool fCopyInput;           // If true, try to copy input file to current directory.
    bool fCopyErrorFatal;      // If true, throw a c++ exception in case of copy errors.

    // Data members.

    art::SourceHelper const& fSourceHelper;
    unsigned int fCurrentEvent;
  };

  typedef art::Source<ubcore::DummySourceDetail> DummySource;

  // Implementations.

  // Constructor.
  DummySourceDetail::DummySourceDetail(fhicl::ParameterSet const& pset,
                                       art::ProductRegistryHelper&,
                                       art::SourceHelper const& helper) :
    fMaxEvents(pset.get<unsigned int>("MaxEvents", 1)),
    fCopyInput(pset.get<bool>("CopyInput", false)),
    fCopyErrorFatal(pset.get<bool>("CopyErrorFatal", true)),
    fSourceHelper(helper),
    fCurrentEvent(0)
  {
    std::cout << "DummySourceDetail constructor called." << std::endl;
    std::cout << "MaxEvents = " << fMaxEvents << std::endl;
    std::cout << "CopyInput = " << fCopyInput << std::endl;
    std::cout << "CopyErrorFatal = " << fCopyErrorFatal << std::endl;
  }

  // Open input file.
  void DummySourceDetail::readFile(std::string const& filename, art::FileBlock*& fb)
  {
    std::cout << "\nDummySourceDetail::readFile called." << std::endl;
    std::cout << "File name: " << filename << std::endl;
    fb = new art::FileBlock(art::FileFormatVersion(1, "DummySource"), filename);

    // Copy input file to current directory?

    if(fCopyInput) {
      std::cout << "Input file will be copied to current directory." << std::endl;

      // Is the input file a url (does it contain a colon)?

      if(filename.find(":") < std::string::npos) {
        std::ostringstream ostr;
        ostr << "Input file is a url: " << filename;
        std::cout << ostr.str() << std::endl;
        if(fCopyErrorFatal)
          throw art::Exception(art::errors::FileOpenError) << ostr.str();
      }
      else {

        // Not a url.

        std::filesystem::path source_path(filename);
        std::string target_path = std::filesystem::current_path() / source_path.filename();
        std::cout << "Target path: " << target_path << std::endl;

        // Check whether source and target are equivalent.

        if(std::filesystem::equivalent(source_path, target_path)) {
          std::cout << "Input file will not be copied because source and target are equivalent."
                    << std::endl;
        }
        else {

          // Do the copy.

          bool copy_ok = false;
          try {
            std::filesystem::copy(source_path, target_path);
            copy_ok = true;
          }
          catch(...) {
            copy_ok = false;
          }
          if(copy_ok)
            std::cout << "Copy successful." << std::endl;
          else {
            std::ostringstream ostr;
            ostr << "Copy failed.";
            std::cout << ostr.str() << std::endl;
            if(fCopyErrorFatal)
              throw art::Exception(art::errors::FileOpenError) << ostr.str();
          }
        }
      }
    }
  }

  // Read event.
  bool DummySourceDetail::readNext(art::RunPrincipal const* const inR,
                                   art::SubRunPrincipal const* const inSR,
                                   art::RunPrincipal*& outR,
                                   art::SubRunPrincipal*& outSR,
                                   art::EventPrincipal*& outE)
  {
    std::cout << "\nDummySourceDetail::readNext called." << std::endl;

    // Return empty event.

    if(fCurrentEvent >= fMaxEvents)
      return false;

    if(inR == 0)
      outR = fSourceHelper.makeRunPrincipal(1, 0);
    if(inSR == 0)
      outSR = fSourceHelper.makeSubRunPrincipal(1, 1, 0);
    outE = fSourceHelper.makeEventPrincipal(1, 1, fCurrentEvent, 0);

    ++fCurrentEvent;

    return true;
  }

  // Close input file.
  void DummySourceDetail::closeCurrentFile()
  {
    std::cout << "\nDummySourceDetail::closeCurrentFile called." << std::endl;
    fCurrentEvent = 0;
  }
}

DEFINE_ART_INPUT_SOURCE(ubcore::DummySource)
