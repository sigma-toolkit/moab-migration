#ifndef MOABREADER_H
#define MOABREADER_H

#include "vtkIOGeometryModule.h"  // For export macro
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkNew.h"  //needed for api signature

class vtkInformation;
class vtkInformationVector;

namespace smoab
{
class Tag;
class Interface;
}  // namespace smoab

class MoabReader : public vtkMultiBlockDataSetAlgorithm
{
  public:
    static MoabReader* New();
    vtkTypeMacro( MoabReader, vtkMultiBlockDataSetAlgorithm ) void PrintSelf( ostream& os, vtkIndent indent );

    // Description:
    // Specify file name of the MOAB mesh file.
    vtkSetStringMacro( FileName );
    vtkGetStringMacro( FileName );

  protected:
    MoabReader();
    ~MoabReader();

    int RequestInformation( vtkInformation* vtkNotUsed( request ),
                            vtkInformationVector** vtkNotUsed( inputVector ),
                            vtkInformationVector* outputVector );

    int RequestData( vtkInformation* vtkNotUsed( request ),
                     vtkInformationVector** vtkNotUsed( inputVector ),
                     vtkInformationVector* outputVector );

  private:
    void CreateSubBlocks( vtkNew< vtkMultiBlockDataSet >& root,
                          smoab::Interface* interface,
                          smoab::Tag const* parentTag,
                          smoab::Tag const* extractTag = NULL );

    void ExtractShell( vtkNew< vtkMultiBlockDataSet >& root, smoab::Interface* interface, smoab::Tag const* parentTag );

    MoabReader( const MoabReader& );   // Not implemented.
    void operator=( const MoabReader& );  // Not implemented.
    char* FileName;
};

#endif  // MoabReader_H
