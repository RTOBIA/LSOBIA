#include "otbWrapperApplication.h"
#include "otbWrapperApplicationFactory.h"

namespace otb
{
namespace Wrapper
{


class LSSegmentation : public otb::Wrapper::Application
{
public:

	typedef LSSegmentation Self;
    typedef Application SuperClass;
	typedef itk::SmartPointer<Self> Pointer;

	itkNewMacro(Self);
	itkTypeMacro(LSSegmentation, otb::Wrapper::Application);

private:
	void DoInit()
	{
		SetName("LSSegmentation");
		SetDescription("Large Scale Image Segmentation Application");

	    //Documentation
	    SetDocName("Large Scale Segmentation");
	    SetDocLongDescription("This application provides several methods to perform segmentation of very high resolution images");
	    SetDocLimitations("None");
	    SetDocAuthors("OBIA-Team");
	    SetDocSeeAlso(" ");

	    // Ram parameter
	    AddRAMParameter();


	    /*
	     * Parameter types :
	     * ParameterType_InputFilenameList "fusion.final.xslist" Final Fusion input XS image list
	     * ParameterType_InputFilename  "fusion.final.pan" "Final Fusion input Panchromatic image "
	     */
	    //
		// add common inputs parameters
		AddParameter(ParameterType_Group,"io","Inputs for Fusion methods (except Final fusion)");
		AddParameter(ParameterType_InputImage,  "io.im",   "Input image");
		SetParameterDescription("io.im", "Image");


	    // output
	    AddParameter(ParameterType_OutputImage, "io.out",  "Result of the segmentation");
	    SetParameterDescription("io.out", "Output Image.");

	    /* CHOICE example
	    AddParameter(ParameterType_Choice,"fusion.final.outformat","Final Fusion output image format");
	    AddChoice("fusion.final.outformat.tiff", "TIFF format" );
	    AddChoice("fusion.final.outformat.jpeg", "JPEG format" );
	    MandatoryOff("fusion.sylvester.linearcombination.image");

	    */

	    /* FLOAT PARAMETER example
	     * AddParameter(ParameterType_Float,"fusion.glp.ratio","Resolutions ratio between the Panchromatic and the multispectral inputs");
    	   SetDefaultParameterFloat("fusion.glp.ratio",  4.);
    	   SetMinimumParameterFloatValue("fusion.bayes.lambda", 0);

	     */


	}

	// Changes the input parameters for some scenarii
	void DoUpdateParameters()
	{
		/*
		 * switch ( GetParameterInt("fusion") )
      {
      case 0:  // finale geometry fusion
      {
		 */
	}

	/*
	 * Ici on met en place le traitement (comme dans un test)
	 */
	void DoExecute()
	{
		int ThisDoesNothing = 0;

	}
};

OTB_APPLICATION_EXPORT(LSSegmentation)
}
}
