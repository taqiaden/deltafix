#include "fixtureModule.hpp"
Expression * fixturemodule::getExpression(char const* nameOfExpression, NXOpen::Part* partName) {
	ExpressionCollection::iterator it1;
	Expression *exp0;
	for (it1 = partName->Expressions()->begin(); it1 != partName->Expressions()->end(); it1++)
	{
		Expression *expression = dynamic_cast<Expression *>(*it1);
		NXString nxname = expression->Name();
		string name = nxname.GetLocaleText();
		if (!strcmp(name.c_str(), nameOfExpression))
		{
			exp0 = expression;
		}
	}
	return exp0;
}
Expression * fixturemodule::getExpression(char const* nameOfExpression) {
	NXOpen::Part *workPart(theSession->Parts()->Work());
	return getExpression(nameOfExpression, workPart);
}
void fixturemodule::messageInfo(char const* message) {
	theUI->NXMessageBox()->Show("", NXOpen::NXMessageBox::DialogTypeError, message);
}
void fixturemodule::messageInfo(string message) {
	theUI->NXMessageBox()->Show("", NXOpen::NXMessageBox::DialogTypeInformation, message.c_str());
}
void fixturemodule::messageInfo(double value) {
	std::ostringstream streamObj;
	streamObj << value;
	theUI->NXMessageBox()->Show("", NXOpen::NXMessageBox::DialogTypeError, (streamObj.str().c_str()));
}
void fixturemodule::messageInfo(NXString value) {
	theUI->NXMessageBox()->Show("", NXOpen::NXMessageBox::DialogTypeError, value.GetText());
}
tag_t fixturemodule::createPoint(double x, double y, double z)
{
	UF_initialize();
	tag_t point_tag = NULL_TAG;
	double point[3] = { x,y,z };//Create a point at x,y,z 
	UF_CSYS_map_point(UF_CSYS_ROOT_WCS_COORDS, point, UF_CSYS_WORK_COORDS, point);//Convert work coordinates to absolute coordinates
	UF_CURVE_create_point(point, &point_tag);//Create point
	UF_terminate();
	return point_tag;
}
tag_t fixturemodule::createPoint(Point3d point3D)
{
	tag_t point_tag = createPoint(point3D.X, point3D.Y, point3D.Z);
	return point_tag;
}
int fixturemodule::checkIfPointinBoundry(double point[3], tag_t body)
{
	int result;
	UF_initialize();
	//1 = point is inside the body
	//2 = point is outside the body
	//3 = point is on the body
	UF_MODL_ask_point_containment(point, body, &result);
	UF_terminate();
	return result;
}
int fixturemodule::checkIfPointinBoundry(Point3d point, tag_t body)
{
	double pointx[3];
	pointx[0] = point.X;
	pointx[1] = point.Y;
	pointx[2] = point.Z;
	return checkIfPointinBoundry(pointx, body);
}
void fixturemodule::highlightObj(tag_t objTag, bool highlight) {
	try{	
		UF_initialize();
		int temp = UF_DISP_set_highlight(objTag, highlight ? 1 : 0);
		UF_DISP_refresh();
		UF_terminate();
	}catch(exception ex)
	{
		UpdateResultSummary("Error: Heighlight attempt failed");
	}
}
void fixturemodule::highlightObj(tag_t objTag) {
	highlightObj(objTag, true);
}
void fixturemodule::DeletePoint(tag_t objTag)
{
	UF_initialize();
	UF_OBJ_delete_object(objTag);
	UF_terminate();
}
void fixturemodule::DeleteAllPoints(vector<tag_t> objTag)
{
	for (int i = 0; i < objTag.size(); i++)
	{
		DeletePoint(objTag[i]);
	}
}
void fixturemodule::drawPoints(int**integerOfPointsCloud, tag_t** &pointsCloud)
{
	for (int i = 0; i < nodesNumber_H; i++) {
		pointsCloud[i] = new tag_t[nodesNumber_V];
		for (int j = 0; j < nodesNumber_V; j++) {
			Point3d currentPoint_transfered;
			Point3d currentPoint_original;
			currentPoint_transfered.X = zeroLocation_transformed.X + (alpha*i);
			currentPoint_transfered.Y= zeroLocation_transformed.Y + (alpha*j);
			currentPoint_transfered.Z= transformed_z_distance;
			currentPoint_original=	ReversedPoint(currentPoint_transfered);
			if (integerOfPointsCloud[i][j] != 2) {
				pointsCloud[i][j] = createPoint(currentPoint_original);
			}
		}
		UF_DISP_refresh();
	}
}
vector<tag_t>  fixturemodule::drawPoints(std::vector<arrayLocation> outputPath)
{
	vector<tag_t> pointsCloud;
	UF_initialize();
	for (int j = 0; j < outputPath.size(); j++) {
		Point3d currentPoint_transfered;
		Point3d currentPoint_original;
		currentPoint_transfered.X= zeroLocation_transformed.X + (alpha*outputPath[j].a);
		currentPoint_transfered.Y = zeroLocation_transformed.Y + (alpha*outputPath[j].b);
		currentPoint_transfered.Z = transformed_z_distance;
		currentPoint_original=	ReversedPoint(currentPoint_transfered);
		pointsCloud.push_back(createPoint(currentPoint_original));
	}
	UF_DISP_refresh();
	UF_terminate();
	return pointsCloud;
}
void fixturemodule::drawPoints(	vector<int> optimumPointsLocation, vector<tag_t> &pointsCloud)
{
	for (int j = 0; j < optimumPointsLocation.size(); j++) {
		vector<	Point3d> locations_transformed(5);
		for (int i = 0; i < locations_transformed.size(); i++)
		{
			locations_transformed[i].X=zeroLocation_transformed.X+(double)selectionDomain[optimumPointsLocation[i]].a*alpha;
			locations_transformed[i].Y=zeroLocation_transformed.Y+(double)selectionDomain[optimumPointsLocation[i]].b*alpha;
			locations_transformed[i].Z=transformed_z_distance;
			Point3d locations_orginal=ReversedPoint(locations_transformed[i]);
			pointsCloud.push_back(createPoint(locations_orginal));
		}
	}
	UF_DISP_refresh();
}
void fixturemodule::refreshPointsHighlight(int**integerOfPointsCloud, tag_t** &pointsCloud)
{
	for (int i = 0; i < nodesNumber_H; i++) {
		for (int j = 0; j < nodesNumber_V; j++) {
			if (integerOfPointsCloud[i][j] != 2 ) {
				highlightObj(pointsCloud[i][j]);
			}
			else {
				highlightObj(pointsCloud[i][j], false);
			}
		}
		UF_MODL_update();
		UF_DISP_refresh();
	}
}
void fixturemodule::refreshPointsHighlight(	vector<arrayLocation>outputPaths, tag_t** &pointsCloud)
{
	for (int i = 0; i < nodesNumber_V; i++) {
		for (int j = 0; j < nodesNumber_H; j++) {
			//Unhighlight All points
			highlightObj(pointsCloud[i][j], false);
		}
	}
	for(int i=0;i<outputPaths.size();i++)
	{
		highlightObj(pointsCloud[outputPaths[i].a][outputPaths[i].b], true);
	}
	UF_MODL_update();
	UF_DISP_refresh();
}
void  fixturemodule::removeAllHighlights(tag_t** pointsCloud)
{
	for (int i = 0; i < nodesNumber_H; i++) {
		for (int j = 0; j < nodesNumber_V; j++) {
			highlightObj(pointsCloud[i][j], false);
		}
	}
}
void fixturemodule::UpdateProgressBar(string info)
{
	UF_initialize();
	progress->SetLabel(info);
	UF_DISP_refresh();
	UF_terminate();
}
void fixturemodule::UpdateResultSummary()
{
	UF_initialize();
	resultSummary->SetValue(result);
	UF_DISP_refresh();
	UF_terminate();
}
void fixturemodule::UpdateResultSummary(string info)
{
	UF_initialize();
	result.push_back(info);
	resultSummary->SetValue(result);
	UF_DISP_refresh();
	UF_terminate();
}
void fixturemodule::UpdateResultSummary_seprator()
{
	UF_initialize();
	result.push_back("-------------------------------------");
	resultSummary->SetValue(result);
	UF_DISP_refresh();
	UF_terminate();
}