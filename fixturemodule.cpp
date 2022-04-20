//==============================================================================
//  WARNING!!  This file is overwritten by the Block UI Styler while generating
//  the automation code. Any modifications to this file will be lost after
//  generating the code again.
//
//      
//
//        This file was generated by the NX Block UI Styler
//        Created by: Taqiaden  AlShameri
//              Version: NX 10
//              Date: 01-18-2021  (Format: mm-dd-yyyy)
//              Time: 14:48 (Format: hh-mm)
//
//==============================================================================
//==============================================================================
//  Purpose:  This TEMPLATE file contains C++ source to guide you in the
//  construction of your Block application dialog. The generation of your
//  dialog file (.dlx extension) is the first step towards dialog construction
//  within NX.  You must now create a NX Open application that
//  utilizes this file (.dlx).
//
//  The information in this file provides you with the following:
//
//  1.  Help on how to load and display your Block UI Styler dialog in NX
//      using APIs provided in NXOpen.BlockStyler namespace
//  2.  The empty callback methods (stubs) associated with your dialog items
//      have also been placed in this file. These empty methods have been
//      created simply to start you along with your coding requirements.
//      The method name, argument list and possible return values have already
//      been provided for you.
//==============================================================================
//------------------------------------------------------------------------------
//These includes are needed for the following template code
//------------------------------------------------------------------------------
#include "fixtureModule.hpp"
#include "calculations.cpp"
using namespace NXOpen;
using namespace NXOpen::BlockStyler;
//------------------------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------------------------
Session *(fixturemodule::theSession) = NULL;
UI *(fixturemodule::theUI) = NULL;
//------------------------------------------------------------------------------
// Constructor for NX Styler class
//------------------------------------------------------------------------------
fixturemodule::fixturemodule()
{
	try
	{
		// Initialize the NX Open C++ API environment
		fixturemodule::theSession = NXOpen::Session::GetSession();
		fixturemodule::theUI = UI::GetUI();
		theDlxFileName = "\\DeltaFix tool files\\fixture module.dlx";
		theDialog = fixturemodule::theUI->CreateDialog(theDlxFileName);
		// Registration of callback functions
		theDialog->AddUpdateHandler(make_callback(this, &fixturemodule::update_cb));
		theDialog->AddCloseHandler(make_callback(this, &fixturemodule::close_cb));
		theDialog->AddFilterHandler(make_callback(this, &fixturemodule::filter_cb));
		theDialog->AddInitializeHandler(make_callback(this, &fixturemodule::initialize_cb));
		theDialog->AddFocusNotifyHandler(make_callback(this, &fixturemodule::focusNotify_cb));
		theDialog->AddKeyboardFocusNotifyHandler(make_callback(this, &fixturemodule::keyboardFocusNotify_cb));
		theDialog->AddDialogShownHandler(make_callback(this, &fixturemodule::dialogShown_cb));
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		throw;
	}
}
//------------------------------------------------------------------------------
// Destructor for NX Styler class
//------------------------------------------------------------------------------
fixturemodule::~fixturemodule()
{
	if (theDialog != NULL)
	{
		delete theDialog;
		theDialog = NULL;
	}
}
//------------------------------- DIALOG LAUNCHING ---------------------------------
//
//    Before invoking this application one needs to open any part/empty part in NX
//    because of the behavior of the blocks.
//
//    Make sure the dlx file is in one of the following locations:
//        1.) From where NX session is launched
//        2.) $UGII_USER_DIR/application
//        3.) For released applications, using UGII_CUSTOM_DIRECTORY_FILE is highly
//            recommended. This variable is set to a full directory path to a file 
//            containing a list of root directories for all custom applications.
//            e.g., UGII_CUSTOM_DIRECTORY_FILE=$UGII_ROOT_DIR\menus\custom_dirs.dat
//
//    You can create the dialog using one of the following way:
//
//    1. USER EXIT
//
//        1) Create the Shared Library -- Refer "Block UI Styler programmer's guide"
//        2) Invoke the Shared Library through File->Execute->NX Open menu.
//
//------------------------------------------------------------------------------
extern "C" DllExport void  ufusr(char *param, int *retcod, int param_len)
{
	fixturemodule *thefixturemodule = NULL;
	try
	{
		thefixturemodule = new fixturemodule();
		// The following method shows the dialog immediately
		thefixturemodule->Show();
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
	if(thefixturemodule != NULL)
	{
		delete thefixturemodule;
		thefixturemodule = NULL;
	}
}
//------------------------------------------------------------------------------
// This method specifies how a shared image is unloaded from memory
// within NX. This method gives you the capability to unload an
// internal NX Open application or user  exit from NX. Specify any
// one of the three constants as a return value to determine the type
// of unload to perform:
//
//
//    Immediately : unload the library as soon as the automation program has completed
//    Explicitly  : unload the library from the "Unload Shared Image" dialog
//    AtTermination : unload the library when the NX session terminates
//
//
// NOTE:  A program which associates NX Open applications with the menubar
// MUST NOT use this option since it will UNLOAD your NX Open application image
// from the menubar.
//------------------------------------------------------------------------------
extern "C" DllExport int ufusr_ask_unload()
{
	//return (int)Session::LibraryUnloadOptionExplicitly;
	return (int)Session::LibraryUnloadOptionImmediately;
	//return (int)Session::LibraryUnloadOptionAtTermination;
}
//------------------------------------------------------------------------------
// Following method cleanup any housekeeping chores that may be needed.
// This method is automatically called by NX.
//------------------------------------------------------------------------------
extern "C" DllExport void ufusr_cleanup(void)
{
	try
	{
		//---- Enter your callback code here -----
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
}
int fixturemodule::Show()
{
	try
	{
		theDialog->Show();
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
	return 0;
}
//------------------------------------------------------------------------------
//---------------------Block UI Styler Callback Functions--------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Callback Name: initialize_cb
//------------------------------------------------------------------------------
void fixturemodule::initialize_cb()
{
	try
	{
		progress = dynamic_cast<NXOpen::BlockStyler::Label*>(theDialog->TopBlock()->FindBlock("progress"));
		group5 = dynamic_cast<NXOpen::BlockStyler::Group*>(theDialog->TopBlock()->FindBlock("group5"));
		drawContactPoints_t = dynamic_cast<NXOpen::BlockStyler::Toggle*>(theDialog->TopBlock()->FindBlock("drawContactPoints_t"));
		showSimulation_t = dynamic_cast<NXOpen::BlockStyler::Toggle*>(theDialog->TopBlock()->FindBlock("showSimulation_t"));
		resultSummary = dynamic_cast<NXOpen::BlockStyler::MultilineString*>(theDialog->TopBlock()->FindBlock("resultSummary"));
		run = dynamic_cast<NXOpen::BlockStyler::Button*>(theDialog->TopBlock()->FindBlock("run"));
		face_select0 = dynamic_cast<NXOpen::BlockStyler::FaceCollector*>(theDialog->TopBlock()->FindBlock("face_select0"));
		group0 = dynamic_cast<NXOpen::BlockStyler::Group*>(theDialog->TopBlock()->FindBlock("group0"));
		iterations_i = dynamic_cast<NXOpen::BlockStyler::IntegerBlock*>(theDialog->TopBlock()->FindBlock("iterations_i"));
		epochs_i = dynamic_cast<NXOpen::BlockStyler::IntegerBlock*>(theDialog->TopBlock()->FindBlock("epochs_i"));
		singulaityDrop = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("singulaityDrop"));
		randomPhasePortion = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("randomPhasePortion"));
		group = dynamic_cast<NXOpen::BlockStyler::Group*>(theDialog->TopBlock()->FindBlock("group"));
		alpha_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("alpha_d"));
		beta_i = dynamic_cast<NXOpen::BlockStyler::IntegerBlock*>(theDialog->TopBlock()->FindBlock("beta_i"));
		group1 = dynamic_cast<NXOpen::BlockStyler::Group*>(theDialog->TopBlock()->FindBlock("group1"));
		useCurvatureCorrection_t = dynamic_cast<NXOpen::BlockStyler::Toggle*>(theDialog->TopBlock()->FindBlock("useCurvatureCorrection_t"));
		zscoreMaxLimit = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("zscoreMaxLimit"));
		increseReactions = dynamic_cast<NXOpen::BlockStyler::Toggle*>(theDialog->TopBlock()->FindBlock("increseReactions"));
		countForFriction_t = dynamic_cast<NXOpen::BlockStyler::Toggle*>(theDialog->TopBlock()->FindBlock("countForFriction_t"));
		frictionCoeficient_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("frictionCoeficient_d"));
		drawingArea0 = dynamic_cast<NXOpen::BlockStyler::DrawingArea*>(theDialog->TopBlock()->FindBlock("drawingArea0"));
		k1_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("k1_d"));
		k2_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("k2_d"));
		k3_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("k3_d"));
		k4_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("k4_d"));
		k5_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("k5_d"));
		k6_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("k6_d"));
		group3 = dynamic_cast<NXOpen::BlockStyler::Group*>(theDialog->TopBlock()->FindBlock("group3"));
		clampingForce_d = dynamic_cast<NXOpen::BlockStyler::DoubleBlock*>(theDialog->TopBlock()->FindBlock("clampingForce_d"));
		machinigForcesData = dynamic_cast<NXOpen::BlockStyler::FileSelection*>(theDialog->TopBlock()->FindBlock("machinigForcesData"));
		ValuesInitilization();
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
}
void fixturemodule::ValuesInitilization()
{
	e=2.71828182845904523536;				//Euler constant
	ebsilon=0.000000000000000001;		//Small value to divide by to avoid singularity when denominator equals zero
	T_min=1; // minimum temperature of DNSA optimizaton algorithm
	T_max=1000;   // maximum temperature of DNSA optimizaton algorithm     
	subGroupSize=10;  
	autoTerminate=false; 	
	sizeOfRandomPhase=300;	
	minimumDiagonalToAlphaRatio=20;
	drawContactPoints=drawContactPoints_t->Value();
	alpha=alpha_d->Value();
	beta=beta_i->Value();
	iterations=iterations_i->Value();
	epochs=epochs_i->Value();
	useCurvetureCorrection=useCurvatureCorrection_t->Value();
	countForFriction=countForFriction_t->Value();
	frictionCoeficient=frictionCoeficient_d->Value();
	clampingForce=clampingForce_d->Value(); 
	k_obj.resize(6);
	k_obj[0]=k1_d->Value();
	k_obj[1]=k2_d->Value();
	k_obj[2]=k3_d->Value();
	k_obj[3]=k4_d->Value();
	k_obj[4]=k5_d->Value();
	k_obj[5]=k6_d->Value();
	showSimulation=showSimulation_t->Value();
	drawContactPoints=drawContactPoints_t->Value();
	meanDistance=100;
	standardDeviaton=1;
	JacobianThreshold=0.0;
	forcesErrorAllowance=0.001;
	frictionOfSmallRegionDiamter=0.3;
	coneAngleInRadian=0.26;    
	numberOfPairs=4;								//Number of adjecent pairs to use when forward and backward interpolating  for approximating first dervative, second dervative and curvature of a point 
	extremPointsAtOuterBoundary=true;  
	dropPercentage=singulaityDrop->Value();
	highZscoreTrigger=zscoreMaxLimit->Value();
	increaseReactionForce=increseReactions->Value();
	maximumBetaToPointsSize=0.2;
	machiningForces_totalExposureScore=0;
	iszeroCuttingForce.clear();
	frictionCoeficient_d->SetEnable(countForFriction_t->Value());
	precentageOfRandomSize=randomPhasePortion->Value();
}
//------------------------------------------------------------------------------
//Callback Name: dialogShown_cb
//This callback is executed just before the dialog launch. Thus any value set 
//here will take precedence and dialog will be launched showing that value. 
//------------------------------------------------------------------------------
void fixturemodule::dialogShown_cb()
{	
	try
	{
		ValuesInitilization();
		UpdateResultSummary();
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
}
//------------------------------------------------------------------------------
//Callback Name: apply_cb
//------------------------------------------------------------------------------
int fixturemodule::apply_cb()
{
	int errorCode = 0;
	try
	{
		//---- Enter your callback code here -----
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		errorCode = 1;
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
	return errorCode;
}
//------------------------------------------------------------------------------
//Callback Name: update_cb
//------------------------------------------------------------------------------
int fixturemodule::update_cb(NXOpen::BlockStyler::UIBlock* block)
{
	try
	{
		if(block == progress)
		{
			//---------Enter your code here-----------
		}
		else if(block == face_select0)
		{
			//---------Enter your code here-----------
			if(			face_select0->GetSelectedObjects().size() >0)
			{
				run->SetEnable(true);
			}
			UF_initialize();
			face=		face_select0->GetSelectedObjects()[0]->Tag();
			UF_SF_face_ask_bounding_box(face,pad_bounding_box);
			UF_terminate();
			//recommend an alpha value
			getTransformedAttributes();
			UpdateResultSummary("Minimum recommended Alpha parameter = "+to_string(boundaryDiagonal/200));
			UpdateResultSummary("Maximum recommended Alpha parameter = "+to_string(boundaryDiagonal/150));
		}
		else if(block == drawContactPoints_t)
		{
			//---------Enter your code here-----------
			drawContactPoints=drawContactPoints_t->Value();
		}
		else if(block == showSimulation_t)
		{
			//---------Enter your code here-----------
			showSimulation=showSimulation_t->Value();
		}
		else if(block == resultSummary)
		{
			//---------Enter your code here-----------
		}
		else if(block == run)
		{
			//---------Enter your code here-----------
			runDNSA();
		}
		else if(block == iterations_i)
		{
			//---------Enter your code here-----------
			iterations=iterations_i->Value();
		}
		else if(block == epochs_i)
		{
			//---------Enter your code here-----------
			epochs=epochs_i->Value();
		}else if(block == singulaityDrop)
		{
			//---------Enter your code here-----------
			dropPercentage=singulaityDrop->Value();
		}
		else if(block == randomPhasePortion)
		{
			//---------Enter your code here-----------
			precentageOfRandomSize=randomPhasePortion->Value();
		}
		else if(block == alpha_d)
		{
			//---------Enter your code here-----------
			alpha=alpha_d->Value();
		}
		else if(block == beta_i)
		{
			//---------Enter your code here-----------
			beta=beta_i->Value();
		}
		else if(block == useCurvatureCorrection_t)
		{
			//---------Enter your code here-----------
			useCurvetureCorrection=useCurvatureCorrection_t->Value();
		} else if(block == zscoreMaxLimit)
		{
			//---------Enter your code here-----------
			highZscoreTrigger=zscoreMaxLimit->Value();
		}
		else if(block == increseReactions)
		{
			//---------Enter your code here-----------
			increaseReactionForce=increseReactions->Value();
		}
		else if(block == countForFriction_t)
		{
			//---------Enter your code here-----------
			countForFriction=countForFriction_t->Value();
			frictionCoeficient_d->SetEnable(countForFriction_t->Value());
		}
		else if(block == frictionCoeficient_d)
		{
			//---------Enter your code here-----------
			frictionCoeficient=frictionCoeficient_d->Value();
		}
		else if(block == drawingArea0)
		{
			//---------Enter your code here-----------
		}
		else if(block == k1_d)
		{
			//---------Enter your code here-----------
			k_obj[0]=k1_d->Value();
		}
		else if(block == k2_d)
		{
			//---------Enter your code here-----------
			k_obj[1]=k2_d->Value();
		}
		else if(block == k3_d)
		{
			//---------Enter your code here-----------
			k_obj[2]=k3_d->Value();
		}
		else if(block == k4_d)
		{
			//---------Enter your code here-----------
			k_obj[3]=k4_d->Value();
		}
		else if(block == k5_d)
		{
			//---------Enter your code here-----------
			k_obj[4]=k5_d->Value();
		}
		else if(block == k6_d)
		{
			//---------Enter your code here-----------
			k_obj[5]=k6_d->Value();
		}
		else if(block == clampingForce_d)
		{
			//---------Enter your code here-----------
			clampingForce=clampingForce_d->Value(); 
		}
		else if(block == machinigForcesData)
		{
			//---------Enter your code here-----------
		}

	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
	return 0;
}
//------------------------------------------------------------------------------
//Callback Name: ok_cb
//------------------------------------------------------------------------------
int fixturemodule::ok_cb()
{
	int errorCode = 0;
	try
	{
		errorCode = apply_cb();
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		errorCode = 1;
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
	return errorCode;
}
//------------------------------------------------------------------------------
//Function Name: GetBlockProperties
//Description: Returns the propertylist of the specified BlockID
//------------------------------------------------------------------------------
PropertyList* fixturemodule::GetBlockProperties(const char *blockID)
{
	return theDialog->GetBlockProperties(blockID);
}
//------------------------------------------------------------------------------
//Callback Name: focusNotify_cb
//This callback is executed when any block (except the ones which receive keyboard entry such as Integer block) receives focus.
//------------------------------------------------------------------------------
void fixturemodule::focusNotify_cb(NXOpen::BlockStyler::UIBlock* block, bool focus)
{
	try
	{
		//---- Enter your callback code here -----
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
}
//------------------------------------------------------------------------------
//Callback Name: keyboardFocusNotify_cb
//This callback is executed when block which can receive keyboard entry, receives the focus.
//------------------------------------------------------------------------------
void fixturemodule::keyboardFocusNotify_cb(NXOpen::BlockStyler::UIBlock* block, bool focus)
{
	try
	{
		//---- Enter your callback code here -----
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
}
//------------------------------------------------------------------------------
//Callback Name: filter_cb
//------------------------------------------------------------------------------
int fixturemodule::filter_cb(NXOpen::BlockStyler::UIBlock* block, NXOpen::TaggedObject* selectObject)
{
	return(UF_UI_SEL_ACCEPT);
}
int fixturemodule::close_cb()
{
	try
	{
		//---- Enter your callback code here -----
	}
	catch(exception& ex)
	{
		//---- Enter your exception handling code here -----
		fixturemodule::theUI->NXMessageBox()->Show("Block Styler", NXOpen::NXMessageBox::DialogTypeError, ex.what());
	}
	return 0;
}
void fixturemodule::clearVariables()
{
	//clear previous records
	result.clear();
	UpdateResultSummary();
	DeleteAllPoints(pointsCloud2_sub); pointsCloud2_sub.clear();
	DeleteAllPoints(pointsCloud3_optimum); pointsCloud3_optimum.clear();
	pointsInForawrdConicalBeam_all.clear();
	lengthOfShearLine_all.clear();
	deflectionScore.clear();
	appliedForcesMoments_atZeroRef.clear();
	pathPointsUnitNorm_y_all.clear();
	pathPointsUnitNorm_x_all.clear();
	real_forcesApplicationPoints_fromZeroReference_Transforemed.clear();
	pointsInForawrdConicalBeam_subSet.clear();
	pointsInForawrdConicalBeam_F.clear();
	pointsInForawrdConicalBeam_all.clear();
	limda1.clear();
	limda2.clear();
	limda3.clear();
	lengthOfShearLine_subSet.clear();
	lengthOfShearLine_F.clear();
	lengthOfShearLine_all.clear();
	deflectionArms.clear();
	areaMomentOfInertia_y.clear();
	areaMomentOfInertia_x.clear();
	selectionDomain.clear();
	pathPointsUnitNorm_y_subSet.clear();
	pathPointsUnitNorm_x_subSet.clear();
	pointsInSmallregion_subSet.clear();
	pointsInBigRegion_subSet.clear();
	pointsInSmallregion_F.clear();
	pointsInBigRegion_F.clear();
	pointsInSmallregion_all.clear();
	pointsInBigRegion_all.clear();
	minimumDistanceToCavity_subSet.clear();
	minimumDistanceToCavity_F.clear();
	minimumDistanceToCavity_all.clear();
	cumulativeDistanceToCavity_subSet.clear();
	cumulativeDistanceToCavity_F.clear();
	cumulativeDistanceToCavity_all.clear();
	JacobianDetermint="";
	detachmentScoreWithFriction="";
	detachmentScoreWithoutFriction="";
	result.clear();
	iszeroCuttingForce.clear();
	forceIsOutsideGeometry.clear();
}
bool fixturemodule::loadANNFiles()
{
	w1=	getCSVFileData("\\DeltaFix tool files\\w1.csv");
	w2=	getCSVFileData("\\DeltaFix tool files\\w2.csv");
	w3=	getCSVFileData("\\DeltaFix tool files\\w3.csv");
	b1=	getCSVFileData("\\DeltaFix tool files\\b1.csv");
	b2=	getCSVFileData("\\DeltaFix tool files\\b2.csv");
	b3=	getCSVFileData("\\DeltaFix tool files\\b3.csv");
	inputMinAndRange=	getCSVFileData("\\DeltaFix tool files\\inputMinAndRange.csv");
	outputMinAndRange=	getCSVFileData("\\DeltaFix tool files\\outputMinAndRange.csv");
	AveDevNormalizationCoeficients=	getCSVFileData("\\DeltaFix tool files\\AveDevNormalizationCoeficients.csv");
	machiningForcesTabular_orginal=	getCSVFileData(machinigForcesData,"\\DeltaFix tool files\\quasi-static loads.csv");
	//Check if files are laoded
	if(w1.size()==0 ||w2.size()==0||w3.size()==0||b1.size()==0||b2.size()==0||b3.size()==0||inputMinAndRange.size()==0||outputMinAndRange.size()==0||AveDevNormalizationCoeficients.size()==0||machiningForcesTabular_orginal.size()==0)
	{
		UpdateResultSummary("Error: Unable to retrieve  the ANN parameters or machining forces tabular data files");
		return false;
	}else
	{
		for (int i = 0; i < machiningForcesTabular_orginal.size(); i++)
		{
			machiningForces_totalExposureScore+=machiningForcesTabular_orginal[i][6];
		}
		for (int i = 0; i < machiningForcesTabular_orginal.size(); i++)
		{
			exposureFriction.push_back(machiningForcesTabular_orginal[i][6]/machiningForces_totalExposureScore);
		}
		return true;
	}
}
std::vector<std::vector<double> >  fixturemodule::getCSVFileData(NXOpen::BlockStyler::FileSelection* blockID, string defaultPath)
{
	std::vector<std::vector<double> > data;
	std::ifstream test(blockID->Path().GetLocaleText()); 
	if( test)
	{
		data=	CSVToArray(blockID->Path().GetLocaleText());	
	}else
	{
		std::ifstream test2(defaultPath); 
		if(test2)
		{
			data=	CSVToArray(defaultPath);
		}
	}
	return data;
}
std::vector<std::vector<double> >  fixturemodule::getCSVFileData(string defaultPath)
{
	std::vector<std::vector<double> > data;
	std::ifstream test2(defaultPath); 
	if(test2)
	{
		data=	CSVToArray(defaultPath);
	}
	return data;
}
void fixturemodule::runDNSA()
{
	//vector<	double> x; x.push_back(9);x.push_back(8);x.push_back(177);x.push_back(1);x.push_back(10);
	//	
	//							
	//							x.erase(x.begin());
	//return;
	//calculation cc;
	//
	//	for (int i = 0; i < 3000; i++)
	//	{
	//		int x=cc.newPointGenerator(0.103906,0.835379,112);
	//		if(x>=112){
	//			messageInfo(	x);
	//		}
	//	}

	//

	////messageInfo(cc.halfRangeShifting(13,13)); //6.5

	//return;



	UpdateProgressBar("Initilize data");
	clearVariables();
	highlightObj(face,false);
	try
	{
		//Record starting time
		auto start = chrono::steady_clock::now();
		calculation cal;
		if(!loadANNFiles())
		{
			return ;
		}
		srand(time(NULL));
		if(!	getTransformedAttributes())
		{
			return;
		}
		if(!	InitilizeForcesMoments_atZeroRef())
		{
			return;
		}
		tag_t** pointsCloud;
		int** integerOfPointsCloud;
		pointsCloud = new tag_t*[nodesNumber_H];
		integerOfPointsCloud = new int*[nodesNumber_H];
		vector<	vector<arrayLocation>> outputPath;
		vector<arrayLocation> extremePoints;
		vector <double>  weights_x;vector <double> weights_y;
		vector <double>curveture_all;
		vector <double>curveture_subSet;
		discretizationUnit( integerOfPointsCloud,outputPath,&extremePoints,&weights_x,&weights_y);
		if(!mathmaticalAttributesUnit(integerOfPointsCloud,outputPath,&extremePoints,&weights_x,&weights_y,curveture_all))
		{
			return;
		}
		UpdateResultSummary("Size of point set domain: "+to_string( curveture_all.size()));
		if(!extractSelectionDomain(curveture_all,curveture_subSet,outputPath))
		{
			return;
		}
		UpdateResultSummary("Size of availabe selection domain: "+to_string( curveture_subSet.size()));
		if(showSimulation)
		{
			UpdateProgressBar("Draw Simulation points");
			pointsCloud2_sub=drawPoints(selectionDomain);
		}
		auto discretizationTime =chrono::steady_clock::now();
		recordElapsedTime( chrono::duration_cast<chrono::seconds>(discretizationTime - start).count());
		string temp=elapsedTime.GetText();
		UpdateResultSummary("Elapsed time = "+ temp);
		vector<int> optimumPointsLocation=		DNSAExecution(&extremePoints,&curveture_subSet,&weights_x,&weights_y,outputPath[0].size());
		auto end =chrono::steady_clock::now();
		recordElapsedTime( chrono::duration_cast<chrono::seconds>(end - start).count());
		if(showSimulation)
		{	
			DeleteAllPoints(pointsCloud2_sub);
		}
		showFinalResult(pointsCloud,optimumPointsLocation);	

		saveStringToFile(& cal.matrixToString(cumulativeDistanceToCavity_subSet),"cumulativeDistanceToCavity_subSet");
		saveStringToFile(& cal.matrixToString(minimumDistanceToCavity_subSet),"minimumDistanceToCavity_subSet");
		saveStringToFile(& cal.matrixToString(pointsInSmallregion_subSet),"pointsInSmallregion_subSet");
		saveStringToFile(& cal.matrixToString(pointsInBigRegion_subSet),"pointsInBigRegion_subSet");
		saveStringToFile(& cal.matrixToString(pathPointsUnitNorm_x_subSet),"pathPointsUnitNorm_x_subSet");
		saveStringToFile(& cal.matrixToString(pathPointsUnitNorm_y_subSet),"pathPointsUnitNorm_y_subSet");
		saveStringToFile(& cal.matrixToString(areaMomentOfInertia_x),"areaMomentOfInertia_x");
		saveStringToFile(& cal.matrixToString(areaMomentOfInertia_y),"areaMomentOfInertia_y");
		saveStringToFile(& cal.matrixToString(lengthOfShearLine_subSet),"lengthOfShearLine_subSet");
		saveStringToFile(& cal.matrixToString(pointsInForawrdConicalBeam_subSet),"pointsInForawrdConicalBeam_subSet");
		saveStringToFile(& cal.matrixToString(curveture_subSet),"curveture_subSet");
		saveStringToFile(& cal.matrixToString(curveture_all),"curveture_all");

		if(drawContactPoints)
		{
			drawPoints(optimumPointsLocation,pointsCloud3_optimum);
		}
		run->SetLabel("Rerun");
	}
	catch(exception& ex)
	{
		UpdateResultSummary( ex.what());
		DeleteAllPoints(pointsCloud2_sub);
		DeleteAllPoints(pointsCloud3_optimum);
	}
}
void fixturemodule::recordElapsedTime(long elapsed_seconds )
{
	double elapsed_seconds_d=elapsed_seconds;
	double elapsed_min=(int)(elapsed_seconds_d/60);
	double elapsed_sec=elapsed_seconds_d-(elapsed_min*60);
	elapsedTime=to_string((int)elapsed_min)+" Min. and "+to_string((int)elapsed_sec)+" sec.";
}