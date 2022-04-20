#include "fixtureModule.hpp"
#include "calculations.cpp"
#include <string>
calculation cal;
Point2d fixturemodule::zeroLocationOfPointsCloud() {
	Point2d result;
	// links is the distance between two nodes
	int numberOfLinks_Horizontal, numberOfLinks_Vertical;
	numberOfLinks_Horizontal = 2 * (int)((transformed_width / 2) / alpha);
	numberOfLinks_Vertical = 2 * (int)((transformed_legnth / 2) / alpha);
	result.X = transformed_CenterPoint.X - ((numberOfLinks_Horizontal / 2)*alpha);
	result.Y = transformed_CenterPoint.Y - ((numberOfLinks_Vertical / 2)*alpha);
	return result;
}
bool fixturemodule::getTransformedAttributes()
{
	UpdateProgressBar("Transform axes");
	//location trasnformation
	transformedAxisMap[0]=0,transformedAxisMap[1]=0,transformedAxisMap[2]=0;
	bool is3dSurface=true;
	for(int i=0;i<3;i++)
	{
		if(abs(pad_bounding_box[i]-pad_bounding_box[i+3])<0.0001)
		{
			transformedAxisMap[2]=i; is3dSurface=false;
		}
	}
	if(is3dSurface)
	{
		messageInfo("The current solution is not available for 3d surfaces");
		return false;
	}
	if(transformedAxisMap[2]==0)
	{
		transformedAxisMap[0]=1;
		transformedAxisMap[1]=2;
	}else	if(transformedAxisMap[2]==1)
	{
		transformedAxisMap[0]=0;
		transformedAxisMap[1]=2;
	}else
	{
		transformedAxisMap[0]=0;
		transformedAxisMap[1]=1;
	}
	reverseAxisMap[0]=0,reverseAxisMap[1]=0,reverseAxisMap[2]=0;
	for(int i=0;i<3;i++)
	{
		reverseAxisMap[transformedAxisMap[i]]=i;
	}
	transformed_z_distance=pad_bounding_box[transformedAxisMap[2]];
	transformed_CenterPoint.X=((pad_bounding_box[transformedAxisMap[0]+3]-pad_bounding_box[transformedAxisMap[0]])/2)+pad_bounding_box[transformedAxisMap[0]];
	transformed_CenterPoint.Y=((pad_bounding_box[transformedAxisMap[1]+3]-pad_bounding_box[transformedAxisMap[1]])/2)+pad_bounding_box[transformedAxisMap[1]];
	//transformed width and legnth
	transformed_width=(pad_bounding_box[transformedAxisMap[0]+3]-pad_bounding_box[transformedAxisMap[0]]);
	transformed_legnth=(pad_bounding_box[transformedAxisMap[1]+3]-pad_bounding_box[transformedAxisMap[1]]);
	nodesNumber_H = (2 * (int)((transformed_width / 2) / alpha)) + 1;UpdateResultSummary("Number of horizontal points : "+to_string(nodesNumber_H));
	nodesNumber_V = (2 * (int)((transformed_legnth / 2) / alpha)) + 1; UpdateResultSummary("Number of vertical points : "+to_string(nodesNumber_V));
	numberOfPointsInBoundedRectangular=(double)nodesNumber_H*(double)nodesNumber_V;
	zeroLocation_transformed = zeroLocationOfPointsCloud();  UpdateResultSummary("Transformed zero location : ("+to_string(zeroLocation_transformed.X)+","+to_string(zeroLocation_transformed.Y)+")");
	boundaryDiagonal=sqrt(transformed_width*transformed_width+transformed_legnth*transformed_legnth);
	if(boundaryDiagonal/alpha<minimumDiagonalToAlphaRatio)
	{
		messageInfo("Please specify a smaller alpha parameter ");
		return false;
	}
	return	true;
}
Point3d fixturemodule::ReversedPoint(Point3d transformedPoint)
{
	double currentPoint_transfered[3];
	currentPoint_transfered[0]=transformedPoint.X;
	currentPoint_transfered[1]=transformedPoint.Y;
	currentPoint_transfered[2]=transformedPoint.Z;
	Point3d currentPoint_original;
	currentPoint_original.X=currentPoint_transfered[reverseAxisMap[0]];
	currentPoint_original.Y=currentPoint_transfered[reverseAxisMap[1]];
	currentPoint_original.Z=currentPoint_transfered[reverseAxisMap[2]];
	return currentPoint_original;
}
Point3d fixturemodule::transformedPoint(Point3d orginalPoint)
{
	double currentPoint_orginal[3];
	currentPoint_orginal[0]=orginalPoint.X;
	currentPoint_orginal[1]=orginalPoint.Y;
	currentPoint_orginal[2]=orginalPoint.Z;
	Point3d currentPoint_transformed;
	currentPoint_transformed.X=currentPoint_orginal[transformedAxisMap[0]];
	currentPoint_transformed.Y=currentPoint_orginal[transformedAxisMap[1]];
	currentPoint_transformed.Z=currentPoint_orginal[transformedAxisMap[2]];
	return currentPoint_transformed;
}
void fixturemodule::faceToPointsCloud(int**integerOfPointsCloud)
{
	numberOfSurfacePoints=0;
	double sum_x=0.0,sum_y=0.0;
	numberOfInCavityPoints=0;
	for (int i = 0; i < nodesNumber_H; i++)
	{
		integerOfPointsCloud[i] = new int[nodesNumber_V];
		bool countCavityTrigger=false;
		int cavityTempCounter=0;
		for (int j = 0; j < nodesNumber_V; j++) 
		{
			if(((i*nodesNumber_V)+j)%50==0)
			{
				double	Progress=((double)(((i*nodesNumber_V)+j+1)/(double)(nodesNumber_V*nodesNumber_H))*100);
				UpdateProgressBar("Face entity discretization progress "+to_string((int)Progress)+"%");
			}
			integerOfPointsCloud[i][j] =0;
			Point3d currentPoint_transfered;
			Point3d currentPoint_original;
			currentPoint_transfered.X = zeroLocation_transformed.X + (alpha*i);
			currentPoint_transfered.Y = zeroLocation_transformed.Y + (alpha*j);
			currentPoint_transfered.Z = transformed_z_distance;
			currentPoint_original=	ReversedPoint(currentPoint_transfered);
			integerOfPointsCloud[i][j] =checkIfPointinBoundry(currentPoint_original, face);		
			if(integerOfPointsCloud[i][j]!=2)
			{
				if(countCavityTrigger==false)
				{
					countCavityTrigger=true;
				}else
				{
					numberOfInCavityPoints+=cavityTempCounter; cavityTempCounter=0;
				}
				sum_x+=i; 	sum_y+=j;
				numberOfSurfacePoints+=1;
			}else
			{
				if(countCavityTrigger)
				{
					cavityTempCounter+=1;
				}
			}
		}
	}
	centroid.X=(sum_x*alpha)/numberOfSurfacePoints; centroid.Y=(sum_y*alpha)/numberOfSurfacePoints;
	//area moment of inertia calculation
	areaMomentOfInertia_x.resize(nodesNumber_H,0);
	areaMomentOfInertia_y.resize(nodesNumber_V,0);
	for (int i = 0; i < nodesNumber_H; i++) {
		for (int j = 0; j < nodesNumber_V; j++) {
			if(integerOfPointsCloud[i][j]!=2)
			{
				double d_aroundX=centroid.Y-j*alpha;
				areaMomentOfInertia_x[i]+=((alpha*alpha*alpha)/12)+(alpha*d_aroundX*d_aroundX); 
				double d_aroundY=centroid.X-i*alpha;
				areaMomentOfInertia_y[j]+=((alpha*alpha*alpha)/12)+(alpha*d_aroundY*d_aroundY);
			}
		}
	}
	for (int i = 0; i < nodesNumber_H; i++) {
		if(areaMomentOfInertia_x[i]==0)
		{
			if(i-1>0)
			{
				areaMomentOfInertia_x[i]=areaMomentOfInertia_x[i-1]/2;
			}else
			{
				areaMomentOfInertia_x[i]=areaMomentOfInertia_x[i+1]/2;
			}
		}
	}
	for (int i = 0; i < nodesNumber_V; i++) {
		if(areaMomentOfInertia_y[i]==0)
		{
			if(i-1>0)
			{
				areaMomentOfInertia_y[i]=areaMomentOfInertia_y[i-1]/2;
			}else
			{
				areaMomentOfInertia_y[i]=areaMomentOfInertia_y[i+1]/2;
			}
		}
	}
}
void fixturemodule::copy(int** &snapShotOfPointsCloud, int** &integerOfPointsCloud) 
{
	snapShotOfPointsCloud = new int*[nodesNumber_H];
	for (int i = 0; i < nodesNumber_H; i++) {
		snapShotOfPointsCloud[i] = new int[nodesNumber_V];
		for (int j = 0; j < nodesNumber_V; j++) {
			snapShotOfPointsCloud[i][j] = integerOfPointsCloud[i][j];
		}
	}
}
int fixturemodule::maximumDistance(vector<arrayLocation> dynamicSelectionSet, vector< arrayLocation> pinLocations)
{
	//loop all selection space and choose the point based on maximum distance summation
	int optimumPointIndex = 0;
	double maxDistanceSum = 0;
	for (int c = 0; c < dynamicSelectionSet.size(); c++) {
		//get the distance sum of the current point
		double tempDistance = 0;
		for (int d = 0; d < pinLocations.size(); d++) {
			tempDistance += cal.distance(pinLocations[d].a, pinLocations[d].b, dynamicSelectionSet[c].a, dynamicSelectionSet[c].b);
		}
		if (tempDistance > maxDistanceSum)
		{
			optimumPointIndex = c;
			maxDistanceSum = tempDistance;
		}
	}
	return optimumPointIndex;
}
bool fixturemodule::checkLocationAllowance(int**PointsCloud) {
	//check if there is still space to select a new location
	for (int i = 0; i < nodesNumber_H; i++) {
		for (int j = 0; j < nodesNumber_V; j++) {
			if (PointsCloud[i][j] != 2) { 	return true;
			}
		}
	}
	return false;
}
void fixturemodule::isolateBorder(int**integerOfPointsCloud)
{
	string er;
	//isolate border
	int nodesNumber_V, nodesNumber_H;
	nodesNumber_H = (2 * (int)((transformed_width / 2) /alpha)) + 1;
	nodesNumber_V = (2 * (int)((transformed_legnth / 2) / alpha)) + 1;
	try{
		bool** vParallel = new bool*[nodesNumber_H];
		for (int i = 0; i < nodesNumber_H; i++) {
			vParallel[i] = new bool[nodesNumber_V];
			for (int j = 0; j < nodesNumber_V; j++) {
				vParallel[i][j] = false;
				er+="i"+to_string(i)+"j"+to_string(j);
				//Increse prefromence by ignoring points outside the regoin
				if (integerOfPointsCloud[i][j] == 2) { continue; }
				//detect inner borders
				if (i + 1 < nodesNumber_H && integerOfPointsCloud[i + 1][j] == 2) {
					vParallel[i][j] = true;
				}
				else if (i - 1 >= 0 && integerOfPointsCloud[i - 1][j] == 2)
				{
					vParallel[i][j] = true;
				}
				else if (j + 1 < nodesNumber_V && integerOfPointsCloud[i][j + 1] == 2)
				{
					vParallel[i][j] = true;
				}
				else if (j - 1 >= 0 && integerOfPointsCloud[i][j - 1] == 2)
				{
					vParallel[i][j] = true;
				}
				//detect outer borders
				if (i + 1 == nodesNumber_H) {
					vParallel[i][j] = true;
				}
				else if (i - 1 < 0)
				{
					vParallel[i][j] = true;
				}
				else if (j + 1 == nodesNumber_V)
				{
					vParallel[i][j] = true;
				}
				else if (j - 1 < 0)
				{
					vParallel[i][j] = true;
				}
			}
		}
		for (int i = 0; i < nodesNumber_H; i++) {
			for (int j = 0; j < nodesNumber_V; j++) {
				if (vParallel[i][j] != true)
				{
					integerOfPointsCloud[i][j] = 2;
				}
			}
		}
	}catch(exception& ex)
	{
		messageInfo(er);
	}
}
void fixturemodule::getBorderPaths(int**integerOfPointsCloud,	vector<	vector<arrayLocation>>&outputPaths)
{
	//The number of closed paths is limited to the value of loopLimit
	for (int loopLimit=0;loopLimit<50;loopLimit++)
	{
		vector<arrayLocation>extractedPath; 
		try{
			bool isOuterPath=false;
			if(loopLimit==0)
			{
				isOuterPath=true;
			}
			extractNewPath(integerOfPointsCloud,extractedPath,isOuterPath);
		}catch(exception& ex)
		{
			break; messageInfo(ex.what());
		}
		if(!extractedPath.empty())
		{	
			outputPaths.push_back(extractedPath);
		}else
		{
			break;
		}
		extractedPath.clear();
		if(!checkLocationAllowance(integerOfPointsCloud))
		{ 
			break;
		}
	}
}
void fixturemodule::extractNewPath(int**integerOfPointsCloud,	vector<arrayLocation>&outputPath,	bool isOuterPath )
{
	bool** vParallel = new bool*[nodesNumber_V];
	int** v2 = new int*[nodesNumber_V];
	arrayLocation pointOne,pointTwo;
	int theta=0;
	for (int i = 0; i < nodesNumber_H; i++) {
		for (int j = 0; j < nodesNumber_V; j++) {
			//find first point in a path
			if (integerOfPointsCloud[i][j]!= 2)
			{ 
				pointOne.a=i; pointOne.b=j;
				outputPath.push_back(pointOne);
				break;
			}  else{ continue; }
		}
		if(!outputPath.empty())
		{
			break;
		}
	}
	//find the whole path
	bool endPathDetector=false;
	int limitGuard=0;
	while (!endPathDetector && limitGuard<70000)
	{
		limitGuard+=1;
		arrayLocation newPointInPath;
		for(int i=0; i<8;i++)
		{
			if(isOuterPath)
			{
				// Actions for outer path
				if(theta==7)
				{
					theta=0;
				}else {
					theta+=1;
				}
				//update  pointTwo for outer path 
				if(theta==0)
				{
					pointTwo.a=pointOne.a-1; pointTwo.b=pointOne.b;
				}else	if(theta==1)
				{
					pointTwo.a=pointOne.a-1; pointTwo.b=pointOne.b-1;
				}else	if(theta==2)
				{
					pointTwo.a=pointOne.a; pointTwo.b=pointOne.b-1;
				}else	if(theta==3)
				{
					pointTwo.a=pointOne.a+1; pointTwo.b=pointOne.b-1;
				}else	if(theta==4)
				{
					pointTwo.a=pointOne.a+1; pointTwo.b=pointOne.b;
				}else	if(theta==5)
				{
					pointTwo.a=pointOne.a+1; pointTwo.b=pointOne.b+1;
				}else	if(theta==6)
				{
					pointTwo.a=pointOne.a; pointTwo.b=pointOne.b+1;
				}else	if(theta==7)
				{
					pointTwo.a=pointOne.a-1; pointTwo.b=pointOne.b+1;
				}
			}else
			{
				// Actions for inner path
				if(theta==0)
				{
					theta=7;
				}else {
					theta-=1;
				}
				//update  pointTwo for inner path
				if(theta==0)
				{
					pointTwo.a=pointOne.a+1; pointTwo.b=pointOne.b;
				}else	if(theta==1)
				{
					pointTwo.a=pointOne.a+1; pointTwo.b=pointOne.b-1;
				}else	if(theta==2)
				{
					pointTwo.a=pointOne.a; pointTwo.b=pointOne.b-1;
				}else	if(theta==3)
				{
					pointTwo.a=pointOne.a-1; pointTwo.b=pointOne.b-1;
				}else	if(theta==4)
				{
					pointTwo.a=pointOne.a-1; pointTwo.b=pointOne.b;
				}else	if(theta==5)
				{
					pointTwo.a=pointOne.a-1; pointTwo.b=pointOne.b+1;
				}else	if(theta==6)
				{
					pointTwo.a=pointOne.a; pointTwo.b=pointOne.b+1;
				}else	if(theta==7)
				{
					pointTwo.a=pointOne.a+1; pointTwo.b=pointOne.b+1;
				}
			}
			if(pointTwo.a>=0 && pointTwo.b>=0 && pointTwo.a<nodesNumber_H && pointTwo.b<nodesNumber_V)
			{
				//check if point is on surface
				if(	integerOfPointsCloud[pointTwo.a][pointTwo.b]!= 2)
				{	
					//procedures when new point in path is found
					newPointInPath=pointTwo;
					//swap pointOne and pointTwo
					arrayLocation tempPoint=pointTwo;
					pointTwo=pointOne;
					pointOne=tempPoint;
					// theta should reflect direction
					theta=(theta+4)%8;
					break;
				}
			}
		}
		//check if the new point is the same as the starting point which implies the end of the path
		if(outputPath[0].a==newPointInPath.a && outputPath[0].b==newPointInPath.b)
		{
			endPathDetector=true;
		}else
		{
			outputPath.push_back(newPointInPath);
		}
	}
	//elmenate all  previous allocated path 
	for(int i=0;i<outputPath.size();i++)
	{
		integerOfPointsCloud[outputPath[i].a][outputPath[i].b]=2;
	}
}
void fixturemodule::getExtremePoints(vector<arrayLocation>&outputPath,vector<arrayLocation>&extremPoints,vector <double> & weights_x,vector <double>& weights_y)
{
	double min_x=outputPath[0].a,max_x=outputPath[0].a,min_y=outputPath[0].b,max_y=outputPath[0].b;
	for (arrayLocation iteratePoints : outputPath)
	{
		if(iteratePoints.a<min_x)
		{
			min_x=iteratePoints.a;
		}
		if(iteratePoints.a>max_x)
		{
			max_x=iteratePoints.a;
		}
		if(iteratePoints.b<min_y)
		{
			min_y=iteratePoints.b;
		}
		if(iteratePoints.b>max_y)
		{
			max_y=iteratePoints.b;
		}
	}
	// point 0 is near min_x,min_y
	//point 1 is near max_x,min_y
	//point 2 is near max_x,max_y
	//point 3 is near  min_x,max_y
	double pointMin_dist[4]={0};
	vector<	arrayLocation> extremePoints_temp;
	for (arrayLocation iteratePoints : outputPath)
	{
		if(extremePoints_temp.empty())
		{
			extremePoints_temp.push_back(iteratePoints);		
			extremePoints_temp.push_back(iteratePoints);
			extremePoints_temp.push_back(iteratePoints);
			extremePoints_temp.push_back(iteratePoints);
			//initialize distance
			pointMin_dist[0]=sqrt((iteratePoints.a-min_x)*(iteratePoints.a-min_x)+(iteratePoints.b-min_y)*(iteratePoints.b-min_y));
			pointMin_dist[1]=sqrt((iteratePoints.a-max_x)*(iteratePoints.a-max_x)+(iteratePoints.b-min_y)*(iteratePoints.b-min_y));
			pointMin_dist[2]=sqrt((iteratePoints.a-max_x)*(iteratePoints.a-max_x)+(iteratePoints.b-max_y)*(iteratePoints.b-max_y));
			pointMin_dist[3]=sqrt((iteratePoints.a-min_x)*(iteratePoints.a-min_x)+(iteratePoints.b-max_y)*(iteratePoints.b-max_y));
			continue;
		}else
		{
			double new_pointMin_dist[4]={0};
			new_pointMin_dist[0]=sqrt((iteratePoints.a-min_x)*(iteratePoints.a-min_x)+(iteratePoints.b-min_y)*(iteratePoints.b-min_y));
			new_pointMin_dist[1]=sqrt((iteratePoints.a-max_x)*(iteratePoints.a-max_x)+(iteratePoints.b-min_y)*(iteratePoints.b-min_y));
			new_pointMin_dist[2]=sqrt((iteratePoints.a-max_x)*(iteratePoints.a-max_x)+(iteratePoints.b-max_y)*(iteratePoints.b-max_y));
			new_pointMin_dist[3]=sqrt((iteratePoints.a-min_x)*(iteratePoints.a-min_x)+(iteratePoints.b-max_y)*(iteratePoints.b-max_y));
			if(new_pointMin_dist[0]<pointMin_dist[0])
			{
				extremePoints_temp[0]=iteratePoints;
				pointMin_dist[0]=new_pointMin_dist[0];
			}
			if(new_pointMin_dist[1]<pointMin_dist[1])
			{
				extremePoints_temp[1]=iteratePoints;
				pointMin_dist[1]=new_pointMin_dist[1];
			}
			if(new_pointMin_dist[2]<pointMin_dist[2])
			{
				extremePoints_temp[2]=iteratePoints;
				pointMin_dist[2]=new_pointMin_dist[2];
			}
			if(new_pointMin_dist[3]<pointMin_dist[3])
			{
				extremePoints_temp[3]=iteratePoints;
				pointMin_dist[3]=new_pointMin_dist[3];
			}
		}
	}
	bool minMaxPointsExsctince[4]={false};
	for(int i=0;i<4;i++)
	{
		if(extremePoints_temp[i].a==min_x)
		{
			minMaxPointsExsctince[0]=true;
		}
		if(extremePoints_temp[i].a==max_x)
		{
			minMaxPointsExsctince[1]=true;
		}
		if(extremePoints_temp[i].b==min_y)
		{
			minMaxPointsExsctince[2]=true;
		}
		if(extremePoints_temp[i].b==max_y)
		{
			minMaxPointsExsctince[3]=true;
		}
	}
	for (arrayLocation iteratePoints : outputPath)
	{
		if(minMaxPointsExsctince[0]==false  &&iteratePoints.a==min_x)
		{
			extremePoints_temp.push_back(iteratePoints);minMaxPointsExsctince[0]=true;
		}
		if(minMaxPointsExsctince[1]==false  &&iteratePoints.a==max_x)
		{
			extremePoints_temp.push_back(iteratePoints);minMaxPointsExsctince[1]=true;
		}
		if(minMaxPointsExsctince[2]==false  &&iteratePoints.b==min_y)
		{
			extremePoints_temp.push_back(iteratePoints);minMaxPointsExsctince[2]=true;
		}
		if(minMaxPointsExsctince[3]==false  &&iteratePoints.b==max_y)
		{
			extremePoints_temp.push_back(iteratePoints);minMaxPointsExsctince[3]=true;
		}
	}
	extremPoints=extremePoints_temp;
	//weights calcuation
	//for x
	vector <double> weights_x_temp(extremePoints_temp.size());
	for (arrayLocation iteratePoints : outputPath)
	{
		vector<	double >nearness(extremePoints_temp.size());
		double sum=0.0;
		for(int i=0;i<extremePoints_temp.size();i++)
		{
			double distance=(double)abs( iteratePoints.a-extremePoints_temp[i].a);
			if(distance==0)
			{
				distance=ebsilon;
			}
			nearness[i]=(double)1/distance;
			sum+=	nearness[i];			
		}
		for(int i=0;i<extremePoints_temp.size();i++)
		{
			weights_x_temp[i]+=	(double)	nearness[i]/sum;
		}
	}
	vector <double> weights_y_temp(extremePoints_temp.size());
	for (arrayLocation iteratePoints : outputPath)
	{
		vector<	double >nearness(extremePoints_temp.size());
		double sum=0.0;
		for(int i=0;i<extremePoints_temp.size();i++)
		{
			double distance=(double)abs( iteratePoints.b-extremePoints_temp[i].b);
			if(distance==0)
			{
				distance=ebsilon;
			}
			nearness[i]=(double)1/distance;
			sum+=	nearness[i];			
		}
		for(int i=0;i<extremePoints_temp.size();i++)
		{
			weights_y_temp[i]+=	(double)	nearness[i]/sum;
		}
	}
	weights_x=weights_x_temp;
	weights_y=weights_y_temp;
}
std::vector<std::vector<double>> fixturemodule::getJacobian(	std::vector <double>*pathPointsUnitNorm_x,	std::vector <double>*pathPointsUnitNorm_y,int locatorsPointsLocation[3],arrayLocation* KPC,	std::vector<arrayLocation>* outputPath)
{
	std::vector<std::vector<double>> J(3);
	for(int i=0;i<3;i++)
	{
		double temp_nx,temp_ny,temp_X,temp_Y;
		temp_nx=(*pathPointsUnitNorm_x)[locatorsPointsLocation[i]];
		temp_ny=(*pathPointsUnitNorm_y)[locatorsPointsLocation[i]];
		temp_X=((*outputPath)[locatorsPointsLocation[i]].a-(*KPC).a)*alpha;
		temp_Y=((*outputPath)[locatorsPointsLocation[i]].b-(*KPC).b)*alpha;
		J[i].push_back(-temp_nx);
		J[i].push_back(-temp_ny);
		J[i].push_back((temp_ny*temp_X)-(temp_nx*temp_Y));
	}
	return J;	
}
bool fixturemodule::InitilizeForcesMoments_atZeroRef()
{
	UpdateProgressBar("Initilize forces and moments");
	bool forceOutsideBoundary=false;
	Point3d temp ;
	Point3d temp_transformed ;
	Point2d temp_transformed_2d_shiftedToZero ;
	for(int i=0;i<machiningForcesTabular_orginal.size();i++)
	{
		temp.X=machiningForcesTabular_orginal[i][3]; temp.Y=machiningForcesTabular_orginal[i][4];	temp.Z=machiningForcesTabular_orginal[i][5];	
		temp_transformed=	transformedPoint(temp);
		temp_transformed_2d_shiftedToZero.X=temp_transformed.X-zeroLocation_transformed.X;
		temp_transformed_2d_shiftedToZero.Y=temp_transformed.Y-zeroLocation_transformed.Y;
		real_forcesApplicationPoints_fromZeroReference_Transforemed.push_back(temp_transformed_2d_shiftedToZero);
	}
	//Check if fources are outside the boundary
	iszeroCuttingForce.resize(machiningForcesTabular_orginal.size(),false);
	for(int i=0;i<machiningForcesTabular_orginal.size();i++)
	{ 
		if( (machiningForcesTabular_orginal[i][0]==0 && machiningForcesTabular_orginal[i][1]==0 && machiningForcesTabular_orginal[i][2]==0))
		{
			iszeroCuttingForce[i]=true;
		}else
		{
			if(real_forcesApplicationPoints_fromZeroReference_Transforemed[i].X<0.0 || 
				real_forcesApplicationPoints_fromZeroReference_Transforemed[i].Y<0.0 ||
				real_forcesApplicationPoints_fromZeroReference_Transforemed[i].X>transformed_width || 
				real_forcesApplicationPoints_fromZeroReference_Transforemed[i].Y>transformed_legnth)
			{
				forceOutsideBoundary=true;
			}
		}
	}
	//define new force moment vector at reference point
	appliedForcesMoments_atZeroRef.resize(machiningForcesTabular_orginal.size(), vector<vector<double>>(3));
	Point3d temp_forces ;
	Point3d temp__forces_transformed ;
	for(int i=0;i<machiningForcesTabular_orginal.size();i++)
	{
		temp_forces.X=machiningForcesTabular_orginal[i][0]; temp_forces.Y=machiningForcesTabular_orginal[i][1]; temp_forces.Z=machiningForcesTabular_orginal[i][2];
		temp__forces_transformed=transformedPoint(temp_forces);
		appliedForcesMoments_atZeroRef[i][0].push_back(temp__forces_transformed.X);
		appliedForcesMoments_atZeroRef[i][1].push_back(temp__forces_transformed.Y);
		appliedForcesMoments_atZeroRef[i][2].push_back(appliedForcesMoments_atZeroRef[i][0][0]*real_forcesApplicationPoints_fromZeroReference_Transforemed[0].Y-
			appliedForcesMoments_atZeroRef[i][1][0]*real_forcesApplicationPoints_fromZeroReference_Transforemed[0].X);
		//View data
		UpdateResultSummary("Force "+to_string(i)+" transformed location is ("+
			to_string(	real_forcesApplicationPoints_fromZeroReference_Transforemed[i].X)+","+
			to_string(	real_forcesApplicationPoints_fromZeroReference_Transforemed[i].Y)+")");
		UpdateResultSummary("Force components, F_x="+to_string(	appliedForcesMoments_atZeroRef[i][0][0])+", F_y="+
			to_string(	appliedForcesMoments_atZeroRef[i][1][0]));
	}
	if(forceOutsideBoundary)
	{
		messageInfo("Force/s are outside the geometry");
		return false;
	}else
	{
		return true;
	}
}
double fixturemodule::costCalculation(	std::vector <double>*pathPointsUnitNorm_x,	std::vector <double>*pathPointsUnitNorm_y,vector<int> clampsAndLocatorsPoints,vector<arrayLocation>* extremePoints,vector<arrayLocation> *selectionDomain,vector <double>*curveture2,	vector <double> * weights_x,vector <double>* weights_y)
{
	string error;
	try{
		int locatorsPointsLocation[3];
		locatorsPointsLocation[0]=clampsAndLocatorsPoints[0];
		locatorsPointsLocation[1]=clampsAndLocatorsPoints[1];
		locatorsPointsLocation[2]=clampsAndLocatorsPoints[2];
		error="0.1";
		int clampsPoint[2];
		clampsPoint[0]=clampsAndLocatorsPoints[3];
		clampsPoint[1]=clampsAndLocatorsPoints[4];
		error="0.2";
		// Normal vectors direction is outward i.e., emerge from the surface
		error="0.2"+to_string( locatorsPointsLocation[0])+","+to_string( locatorsPointsLocation[1])+","+to_string( locatorsPointsLocation[2]);
		vector<vector	<	double>> Nt=	cal.getNt(pathPointsUnitNorm_x,pathPointsUnitNorm_y,locatorsPointsLocation);
		error="0.3";
		vector<vector<vector<double>>> J((*extremePoints).size(), vector<vector<double>>(3));
		error="0.4";
		for(int z=0;z<J.size();z++)
		{
			J[z]=	getJacobian(pathPointsUnitNorm_x,pathPointsUnitNorm_y,locatorsPointsLocation,&(*extremePoints)[z],selectionDomain);
		}
		error="1";
		int matrixRank;
		if(abs( cal.determinant(J[0]))<=JacobianThreshold)
		{
			localizationCost=0;
			return 0;
		}
		JacobianDetermint+=to_string( cal.determinant(J[0]))+"\n";
		matrixRank=cal.compute_rank(J[0]);	
		if(matrixRank!=3)
		{
			localizationCost=0;
			return 0;
		}
		error="2";
		////compute reactions
		//Jacobian at reference (0,0)
		vector<vector<double>>	forceJacobian(3);
		arrayLocation forceApplicationPoints={0,0}; 
		forceJacobian=getJacobian(pathPointsUnitNorm_x,pathPointsUnitNorm_y,locatorsPointsLocation,&forceApplicationPoints,selectionDomain);
		error="3";
		//define resultant forces
		vector<vector<vector<double>>> reultantCuttingAndClampingForces(appliedForcesMoments_atZeroRef.size(), vector<vector<double>>(3));
		double hypotenuse_clamp1=sqrt(pathPointsUnitNorm_x[0][clampsPoint[0]]*pathPointsUnitNorm_x[0][clampsPoint[0]]+pathPointsUnitNorm_y[0][clampsPoint[0]]*pathPointsUnitNorm_y[0][clampsPoint[0]]);
		double hypotenuse_clamp2=sqrt(pathPointsUnitNorm_x[0][clampsPoint[1]]*pathPointsUnitNorm_x[0][clampsPoint[1]]+pathPointsUnitNorm_y[0][clampsPoint[1]]*pathPointsUnitNorm_y[0][clampsPoint[1]]);
		// we add the minius because the norms are emerged outward from the workpiece
		double clamp1Force_x=(-(*pathPointsUnitNorm_x)[clampsPoint[0]]*clampingForce)/hypotenuse_clamp1;
		double clamp2Force_x=(-(*pathPointsUnitNorm_x)[clampsPoint[1]]*clampingForce)/hypotenuse_clamp2;
		double clamp1Force_y=(-(*pathPointsUnitNorm_y)[clampsPoint[0]]*clampingForce)/hypotenuse_clamp1;
		double clamp2Force_y=(-(*pathPointsUnitNorm_y)[clampsPoint[1]]*clampingForce)/hypotenuse_clamp2;
		double clampingTotalForce_x=clamp1Force_x+clamp2Force_x;
		double clampingTotalForce_y=clamp1Force_y+clamp2Force_y;
		double clamp1MomentAtZeroRef=clamp1Force_x*(*selectionDomain)[clampsPoint[0]].b*alpha-clamp1Force_y*(*selectionDomain)[clampsPoint[0]].a*alpha;
		double clamp2MomentAtZeroRef=clamp2Force_x*(*selectionDomain)[clampsPoint[1]].b*alpha-clamp2Force_y*(*selectionDomain)[clampsPoint[1]].a*alpha;
		double clampsTotalMomentAtZeroRef=clamp1MomentAtZeroRef+clamp2MomentAtZeroRef;
		error="4";
		// resultant cutting and clamping forces at zero reference point
		for(int i=0;i<appliedForcesMoments_atZeroRef.size();i++)
		{
			reultantCuttingAndClampingForces[i][0].push_back(appliedForcesMoments_atZeroRef[i][0][0]+clampingTotalForce_x);
			reultantCuttingAndClampingForces[i][1].push_back(appliedForcesMoments_atZeroRef[i][1][0]+clampingTotalForce_y);
			reultantCuttingAndClampingForces[i][2].push_back(appliedForcesMoments_atZeroRef[i][2][0]+clampsTotalMomentAtZeroRef);
		}
		error="5";
		//Reaction force matrix
		vector<vector<vector<double>>> locatorsNormalForces(reultantCuttingAndClampingForces.size(), vector<vector<double>>(3));
		vector<vector<vector<double>>> locatorsFrictionForces(reultantCuttingAndClampingForces.size(), vector<vector<double>>(3));
		vector<vector<vector<double>>> limda(reultantCuttingAndClampingForces.size());
		vector<vector<vector<double>>> limda_withoutFriction(reultantCuttingAndClampingForces.size());  
		error="6";
		for(int i=0;i<reultantCuttingAndClampingForces.size();i++)
		{
			locatorsFrictionForces[i][0].push_back(0); // x friction at locator 1
			locatorsFrictionForces[i][0].push_back(0); // y friction at locator 1
			locatorsFrictionForces[i][1].push_back(0); // x friction at locator 2
			locatorsFrictionForces[i][1].push_back(0); // y friction at locator 2
			locatorsFrictionForces[i][2].push_back(0); // x friction at locator 3
			locatorsFrictionForces[i][2].push_back(0); // y friction at locator 3
			error="6.1";
			//Find Normal forces
			//zhama matrix
			vector<vector<double>> Jacobian_inve_trans_zeroRef=cal.getTranspose( cal.getInverse(forceJacobian)); 
			error="6.1.1";
			limda[i]=cal.multiply(Jacobian_inve_trans_zeroRef,reultantCuttingAndClampingForces[i]); 
			error="6.1.2";
			limda_withoutFriction[i]=limda[i]; //keep data for recording purpose 
			error="6.2";
			locatorsNormalForces[i][0].push_back(Nt[0][0]*limda[i][0][0]); //x norml force at locator 1
			locatorsNormalForces[i][0].push_back(Nt[0][1]*limda[i][0][0]); //y norml force at locator 1
			locatorsNormalForces[i][1].push_back(Nt[1][2]*limda[i][1][0]);
			locatorsNormalForces[i][1].push_back(Nt[1][3]*limda[i][1][0]);
			locatorsNormalForces[i][2].push_back(Nt[2][4]*limda[i][2][0]);
			locatorsNormalForces[i][2].push_back(Nt[2][5]*limda[i][2][0]);
			error="6.3";
			if(countForFriction)
			{
				vector<vector<double>> locatorsNormalForces_modified=locatorsNormalForces[i];
				double cumulativeError;
				for(int s=0;s<300;s++)
				{
					//friction direction
					//estimate new values for friction
					locatorsFrictionForces[i][0][0]=-(frictionCoeficient*locatorsNormalForces_modified[0][1]); // x friction at locator 1
					locatorsFrictionForces[i][0][1]=-(frictionCoeficient*locatorsNormalForces_modified[0][0]); // y friction at locator 1
					locatorsFrictionForces[i][1][0]=-(frictionCoeficient*locatorsNormalForces_modified[1][1]); // x friction at locator 2
					locatorsFrictionForces[i][1][1]=-(frictionCoeficient*locatorsNormalForces_modified[1][0]); // y friction at locator 2
					locatorsFrictionForces[i][2][0]=-(frictionCoeficient*locatorsNormalForces_modified[2][1]); // x friction at locator 3
					locatorsFrictionForces[i][2][1]=-(frictionCoeficient*locatorsNormalForces_modified[2][0]); // y friction at locator 3
					vector<vector<double>> modifiedForcesMoments(3);
					modifiedForcesMoments[0].push_back(reultantCuttingAndClampingForces[i][0][0]+locatorsFrictionForces[i][0][0]+locatorsFrictionForces[i][1][0]+locatorsFrictionForces[i][2][0]);
					modifiedForcesMoments[1].push_back(reultantCuttingAndClampingForces[i][1][0]+locatorsFrictionForces[i][0][1]+locatorsFrictionForces[i][1][1]+locatorsFrictionForces[i][2][1]);
					//moment is relative to zer zero reference point
					double xForce_yArm=	locatorsFrictionForces[i][0][0]*	(*selectionDomain)[locatorsPointsLocation[0]].b*alpha	+
						locatorsFrictionForces[i][1][0]*	(*selectionDomain)[locatorsPointsLocation[1]].b*alpha	+
						locatorsFrictionForces[i][2][0]*	(*selectionDomain)[locatorsPointsLocation[2]].b*alpha	;
					double yForce_xArm=	locatorsFrictionForces[i][0][1]*	(*selectionDomain)[locatorsPointsLocation[0]].a*alpha	+
						locatorsFrictionForces[i][1][1]*	(*selectionDomain)[locatorsPointsLocation[1]].a*alpha	+
						locatorsFrictionForces[i][2][1]*	(*selectionDomain)[locatorsPointsLocation[2]].a*alpha	;
					modifiedForcesMoments[2].push_back(reultantCuttingAndClampingForces[i][2][0]+xForce_yArm-yForce_xArm);
					vector<vector<double>> limda_modified=cal.multiply(Jacobian_inve_trans_zeroRef,modifiedForcesMoments);
					//update cumulative error
					cumulativeError=0;
					if(	locatorsNormalForces_modified[0][0]!=0)
					{
						cumulativeError+=abs((locatorsNormalForces_modified[0][0]-(Nt[0][0]*limda_modified[0][0]))/locatorsNormalForces_modified[0][0]);
					}
					if(	locatorsNormalForces_modified[0][1]!=0)
					{
						cumulativeError+=abs((locatorsNormalForces_modified[0][1]-(Nt[0][1]*limda_modified[0][0]))/locatorsNormalForces_modified[0][1]);
					}
					if(	locatorsNormalForces_modified[1][0]!=0)
					{
						cumulativeError+=abs((locatorsNormalForces_modified[1][0]-(Nt[1][2]*limda_modified[1][0]))/locatorsNormalForces_modified[1][0]);
					}
					if(	locatorsNormalForces_modified[1][1]!=0)
					{
						cumulativeError+=abs((locatorsNormalForces_modified[1][1]-(Nt[1][3]*limda_modified[1][0]))/locatorsNormalForces_modified[1][1]);
					}
					if(	locatorsNormalForces_modified[2][0]!=0)
					{
						cumulativeError+=abs((locatorsNormalForces_modified[2][0]-(Nt[2][4]*limda_modified[2][0]))/locatorsNormalForces_modified[2][0]);
					}
					if(	locatorsNormalForces_modified[2][1]!=0)
					{
						cumulativeError+=abs((locatorsNormalForces_modified[2][1]-(Nt[2][5]*limda_modified[2][0]))/locatorsNormalForces_modified[2][1]);
					}
					if(s!=0 && cumulativeError<forcesErrorAllowance)
					{ 
						locatorsNormalForces[i]=locatorsNormalForces_modified;
						limda[i]=limda_modified;
						break;
					}
					locatorsNormalForces_modified[0][0]=(locatorsNormalForces_modified[0][0]+(Nt[0][0]*limda_modified[0][0]))/2; //x norml force at locator 1
					locatorsNormalForces_modified[0][1]=(locatorsNormalForces_modified[0][1]+(Nt[0][1]*limda_modified[0][0]))/2; //y norml force at locator 1
					locatorsNormalForces_modified[1][0]=(locatorsNormalForces_modified[1][0]+(Nt[1][2]*limda_modified[1][0]))/2;
					locatorsNormalForces_modified[1][1]=(locatorsNormalForces_modified[1][1]+(Nt[1][3]*limda_modified[1][0]))/2;
					locatorsNormalForces_modified[2][0]=(locatorsNormalForces_modified[2][0]+(Nt[2][4]*limda_modified[2][0]))/2;
					locatorsNormalForces_modified[2][1]=(locatorsNormalForces_modified[2][1]+(Nt[2][5]*limda_modified[2][0]))/2;
				}
			}
		}
		error="7";
		detachmentScore=0.0;
		locatorsReactions=0.0;
		limdas_spreed=0.0;
		limda1.clear();limda2.clear();limda3.clear();
		//constants for detachment penalty
		double a=(3*clampingForce)/10;
		double ee=((45*a*a)-(9*a))/((41*a)+8);

		for(int z=0;z<limda.size();z++)
		{
			for (int i = 0; i < 3; i++)
			{
				double reversedLimda=-limda[z][i][0];
				if(reversedLimda<=0)
				{
					detachmentScore+=exposureFriction[z]*(abs(reversedLimda)+(10*a));
				}else{ 
					if(increaseReactionForce)
					{
						locatorsReactions+=exposureFriction[z]*(1/reversedLimda);
					}else
					{
						locatorsReactions+=exposureFriction[z]*reversedLimda;
					}

					if(reversedLimda<=a)
					{
						detachmentScore+=exposureFriction[z]*((10*a)-reversedLimda-ee);
					}else if(reversedLimda<=5*a)
					{
						detachmentScore+=exposureFriction[z]*(( ((45*a*a)-(9*a)) / ((9*reversedLimda)-(4*a)+8) )-ee);
					}
				} 
			}
			//UpdateResultSummary("detachment for z"+to_string(z)+" = "+to_string(detachmentScore)+" exposure = "+to_string(exposureFriction[z]));
			limda1.push_back(limda[z][0][0]); limda2.push_back(limda[z][1][0]);  limda3.push_back(limda[z][2][0]); 
			vector<double> currentLimdas; currentLimdas.push_back(limda[z][0][0]); currentLimdas.push_back(limda[z][1][0]); currentLimdas.push_back(limda[z][2][0]);
			limdas_spreed+=exposureFriction[z]*(pow(abs(currentLimdas[0]-currentLimdas[1]),2)+pow(abs(currentLimdas[0]-currentLimdas[2]),2)+pow(abs(currentLimdas[1]-currentLimdas[2]),2));
		}
		if(minimumrDetachmentScore==0 || minimumrDetachmentScore>detachmentScore)
		{
			minimumrDetachmentScore=detachmentScore;
		}
		error="8";
		if(countForFriction)
		{
			detachmentScoreWithFriction+=to_string(detachmentScore)+"\n";
		}
		double detachmentScore_withoutFriction=0.0;
		for(int z=0;z<limda_withoutFriction.size();z++)
		{
			if(limda_withoutFriction[z][0][0]>0)
			{
				detachmentScore_withoutFriction+=limda_withoutFriction[z][0][0];
			}
			if(limda_withoutFriction[z][1][0]>0)
			{detachmentScore_withoutFriction+=limda_withoutFriction[z][1][0];
			}
			if(limda_withoutFriction[z][2][0]>0)
			{
				detachmentScore_withoutFriction+=limda_withoutFriction[z][2][0];
			}
		}
		detachmentScoreWithoutFriction+=to_string(detachmentScore_withoutFriction)+"\n";
		error="9";
		vector<		vector<vector	<	double>>> JInverse((*extremePoints).size(), vector<vector<double>>(3));
		for(int z=0;z<J.size();z++)
		{
			if(matrixRank==3)
			{
				JInverse[z]=		cal.getInverse(J[z]);
			}
		}
		error="10";
		vector<	double> x_dev(J.size()),y_dev(J.size()),theta_dev(J.size());
		for(int z=0;z<J.size();z++)
		{
			vector<vector	<	double>> temp_neg_JInverse_Nt;
			temp_neg_JInverse_Nt=cal.negative(JInverse[z]);
			temp_neg_JInverse_Nt=cal.multiply(temp_neg_JInverse_Nt,Nt);
			//shortcut
			double varience=standardDeviaton*standardDeviaton;
			vector<vector	<	double>> R(6);
			for(int i=0;i<R.size();i++)
			{
				int locator=locatorsPointsLocation[(int)(0.5+(2*((double)i/5)))];
				if(useCurvetureCorrection)
				{
					R[i].push_back(varience*((1+(meanDistance*abs((*curveture2)[locator])))) );
				}else
				{
					R[i].push_back(varience );
				}
			}
			error="11";
			cal.squareAllMatrixMembers(temp_neg_JInverse_Nt);
			std::vector<std::vector<double>>  varienceCoeficient=		cal.multiply(	temp_neg_JInverse_Nt,R);
			x_dev[z]=(*weights_x)[z]* varienceCoeficient[0][0];
			y_dev[z]=(*weights_y)[z]*varienceCoeficient[1][0];
		}
		localizationCost=0;
		for(int z=0;z<J.size();z++)
		{
			localizationCost+=x_dev[z];
			localizationCost+=y_dev[z];
		}
		error="12";
		if(minimumLocaliztionScore==0 || minimumLocaliztionScore>localizationCost)
		{
			minimumLocaliztionScore=localizationCost;
		}
		frictionOfLoclizationObject=0; frictionOfrestrainObject=1;
		double totalCost= frictionOfLoclizationObject*(localizationCost/1) + frictionOfrestrainObject*(detachmentScore/1);
		return (localizationCost*(1+ detachmentScore));
	}catch(exception ex)
	{
		result.push_back("evaluation error level is at "+error);
		UpdateResultSummary();
	}
}
bool fixturemodule::Dervatives(vector<	vector<arrayLocation>>* outputPath,int numberOfPairs,double ebsilon,vector <Point2d>&backpoint,vector <Point2d>&forwardpoint,vector <double>&firstDervative2,	vector <double>&secondDervative2,vector <double>&curveture2)
{
	string errorLevel;
	try{
		//Find the average distance between all pairs of two adjusent points
		double distance_sum=0;
		for(int i=0;i<(*outputPath)[0].size();i++)
		{
			errorLevel="1";
			int previousPoint;
			if(i==0)
			{
				previousPoint=(*outputPath)[0].size()-1;
			}else
			{
				previousPoint=i-1;
			}
			errorLevel="2";
			double x_dist=(*outputPath)[0][i].a- (*outputPath)[0][previousPoint].a ;
			double y_dist=(*outputPath)[0][i].b- (*outputPath)[0][previousPoint].b ;
			distance_sum+=alpha* sqrt((x_dist*x_dist)+(y_dist*y_dist) );
			errorLevel="3";
		}
		double d=distance_sum/(*outputPath)[0].size(); 
		errorLevel="4";
		for(int i=0;i<(*outputPath)[0].size();i++)
		{
			if(i%10==0)
			{
				double	Progress=((double)((i)/(double)((*outputPath)[0].size()))*100);
				UpdateProgressBar("Compute dervatives progress "+to_string((int)Progress)+"%");
			}
			errorLevel="5";
			std::vector<double> backwardPoints_x,backwardPoints_y;
			std::vector<double> forwardPoints_x,forwardPoints_y;
			std::vector<double> back_Points_x,back_Points_y;
			for(int j=0;j<numberOfPairs;j++)
			{
				errorLevel="6";
				//next value
				forwardPoints_x.push_back((*outputPath)[0][(i+j)%(*outputPath)[0].size()].a);
				forwardPoints_y.push_back((*outputPath)[0][(i+j)%(*outputPath)[0].size()].b);
				//previous value
				int position_temp=i-j;
				errorLevel="7";
				if(position_temp>=0)
				{
					errorLevel="8";
					backwardPoints_x.push_back((*outputPath)[0][position_temp].a);
					backwardPoints_y.push_back((*outputPath)[0][position_temp].b);
				}else
				{	
					errorLevel="9";
					backwardPoints_x.push_back((*outputPath)[0][(*outputPath)[0].size()+position_temp].a);
					backwardPoints_y.push_back((*outputPath)[0][(*outputPath)[0].size()+position_temp].b);
				}
			}
			errorLevel="10";
			//backPoint
			Point2d temp;
			double xDif_sum=0,yDif_sum=0,m=0;
			for(int z=1;z<backwardPoints_x.size();z++)
			{	
				errorLevel="11";
				yDif_sum+=(backwardPoints_x.size()-z)*(backwardPoints_y[z]-backwardPoints_y[0]);
				xDif_sum+=(backwardPoints_x.size()-z)*(backwardPoints_x[z]-backwardPoints_x[0]);
			}
			errorLevel="12";
			if(xDif_sum==0)
			{
				xDif_sum=ebsilon;
			}
			m=yDif_sum/xDif_sum;
			temp.X=sqrt((d*d)/((m*m)+1));
			temp.Y=m*temp.X;
			errorLevel="3";
			if(xDif_sum>=0)
			{
				errorLevel="14";
				temp.X=abs(temp.X);
			}else{
				temp.X=-abs(temp.X);
			}
			if(yDif_sum>=0)
			{
				temp.Y=abs(temp.Y);
			}else
			{
				temp.Y=-abs(temp.Y);
			}
			errorLevel="15";
			temp.X+=backwardPoints_x[0];
			temp.Y+=backwardPoints_y[0];
			backpoint.push_back(temp);
			errorLevel="16";
			//forward point
			temp.X=0;			 temp.Y=0;
			xDif_sum=0,yDif_sum=0,m=0;
			for(int z=1;z<forwardPoints_x.size();z++)
			{
				errorLevel="17";
				yDif_sum+=(forwardPoints_x.size()-z)*(forwardPoints_y[z]-forwardPoints_y[0]);
				xDif_sum+=(forwardPoints_x.size()-z)*(forwardPoints_x[z]-forwardPoints_x[0]);
			}
			errorLevel="18";
			if(xDif_sum==0)
			{
				xDif_sum=ebsilon;
			}
			errorLevel="19";
			m=yDif_sum/xDif_sum;
			temp.X=sqrt((d*d)/((m*m)+1));
			temp.Y=m*temp.X;
			errorLevel="20";
			if(xDif_sum>=0)
			{
				errorLevel="21";
				temp.X=abs(temp.X);
			}else{
				temp.X=-abs(temp.X);
			}
			if(yDif_sum>=0)
			{
				errorLevel="22";
				temp.Y=abs(temp.Y);
			}else{
				temp.Y=-abs(temp.Y);
			}
			temp.X+=forwardPoints_x[0];
			temp.Y+=forwardPoints_y[0];
			forwardpoint.push_back(temp);
			double dy,dx;
			dy=forwardpoint[i].Y-backpoint[i].Y;
			dx=forwardpoint[i].X-backpoint[i].X;
			errorLevel="23";
			if(dx==0)
			{
				dx=ebsilon;
			}
			firstDervative2.push_back(dy/dx);
			dy=0;dx=0;
			dy=forwardpoint[i].Y-(2*(*outputPath)[0][i].b) +backpoint[i].Y;
			dx=(forwardpoint[i].X-backpoint[i].X)*(forwardpoint[i].X-backpoint[i].X);
			errorLevel="24";
			if(dx==0)
			{
				dx=ebsilon;
			}
			secondDervative2.push_back(dy/dx);
			dy=0;dx=0;
			double signCheck=(((*outputPath)[0][i].a-backpoint[i].X)*(forwardpoint[i].Y-backpoint[i].Y))		-		(((*outputPath)[0][i].b-backpoint[i].Y)*(forwardpoint[i].X-backpoint[i].X));
			errorLevel="25";
			if(signCheck>=0)
			{
				dy=abs(secondDervative2[i]);
			}else
			{
				dy=-abs(secondDervative2[i]);
			}
			errorLevel="26";
			dx=pow( 1+(firstDervative2[i]*firstDervative2[i]),3/2);
			if(dx==0)
			{
				dx=ebsilon;
			}
			errorLevel="27";
			curveture2.push_back(dy/dx);
		}
		return true;
	}catch(exception ex)
	{
		UpdateResultSummary("Error when computing dervatives");
		UpdateResultSummary("Error level ->"+errorLevel);
		return false;
	}
}
std::vector<std::vector<double> > fixturemodule::CSVToArray(string CSVpath)
{
	std::ifstream  data(CSVpath);
	std::string line;
	std::vector<std::vector<double> > parsedCsv;
	char *end;
	while(std::getline(data,line))
	{
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<double> parsedRow;
		while(std::getline(lineStream,cell,','))
		{
			parsedRow.push_back(strtod(cell.c_str(), &end));
		}
		parsedCsv.push_back(parsedRow);
	}
	data.close();
	return	parsedCsv;
};
double fixturemodule::annOutput(vector< vector<double>>* input)
{
	std::vector<std::vector<double> > Input_AveDevNormalization(1);
	for(int i=0;i<(*input)[0].size();i++)
	{
		//at i==0 is held for output
		Input_AveDevNormalization[0].push_back(((*input)[0][i]-AveDevNormalizationCoeficients[0][i+1])/AveDevNormalizationCoeficients[1][i+1]);
	}
	std::vector<std::vector<double> > Input_MinMaxNormalization(1);
	for(int i=0;i<Input_AveDevNormalization[0].size();i++)
	{
		Input_MinMaxNormalization[0].push_back(((Input_AveDevNormalization[0][i]-inputMinAndRange[0][i])/inputMinAndRange[1][i])*2-1);
	}
	std::vector<std::vector<double> > w1_x(w1.size(),std::vector<double>(w1[0].size()));
	for(int i=0;i<w1.size();i++)
	{
		for(int j=0;j<w1[0].size();j++)
		{
			w1_x[i][j]=w1[i][j]*Input_MinMaxNormalization[0][j];
		}
	}
	std::vector<std::vector<double> > w2_x(w1.size());
	for (int i = 0; i < w1.size(); i++)
	{
		double wieghtedXAndBias_sum=0.0;
		for(int j=0;j<w1[0].size();j++)
		{
			wieghtedXAndBias_sum+=w1_x[i][j];
		}
		wieghtedXAndBias_sum+=b1[i][0];
		double tansig=(2/(1+pow(e,-2*wieghtedXAndBias_sum)))-1;
		w2_x[i].push_back(tansig*w2[0][i]);
	}
	double normalizedOutput=0.0;
	for (int i = 0; i < w2_x.size(); i++)
	{
		normalizedOutput+=w2_x[i][0];
	}
	normalizedOutput+=b2[0][0];
	//reverse min max normalization
	double outputNormlized=(((normalizedOutput+1)/2)*outputMinAndRange[1][0])+outputMinAndRange[0][0];
	//reverse ave dev normalizaton
	double	output=(outputNormlized*AveDevNormalizationCoeficients[1][0])+AveDevNormalizationCoeficients[0][0];
	return output;
}
double fixturemodule::annOutput2(vector< vector<double>>* input)
{
	//network specifications
	/*2 hidden layers	
	transfer function 1	radbasn
	transfer function 2	radbasn
	transfer function 3	purlin
	# neurons hidden layer 1	20
	# neurons hidden layer 2	5*/
	//matlab building code	
	/*tf='radbasn';
	net2 = feedforwardnet([20 5]);
	net2.trainFcn = 'trainbfg';
	net2.layers{1}.transferFcn=tf;
	net2.layers{2}.transferFcn=tf;
	net2=train(net2,in,out);
	y=sim(net2,test);*/
	//z-score standarization
	std::vector<std::vector<double> > Input_AveDevNormalization(1);
	for(int i=0;i<(*input)[0].size();i++)
	{
		//at i==0 is held for output
		Input_AveDevNormalization[0].push_back(((*input)[0][i]-AveDevNormalizationCoeficients[0][i+1])/AveDevNormalizationCoeficients[1][i+1]);
	}
	//min max normalization
	std::vector<std::vector<double> > Input_MinMaxNormalization(1);
	for(int i=0;i<Input_AveDevNormalization[0].size();i++)
	{
		Input_MinMaxNormalization[0].push_back(((Input_AveDevNormalization[0][i]-inputMinAndRange[0][i])/inputMinAndRange[1][i])*2-1);
	}
	//w1*x1
	std::vector<std::vector<double> > w1_x(w1.size(),std::vector<double>(w1[0].size()));
	for(int i=0;i<w1.size();i++)
	{
		for(int j=0;j<w1[0].size();j++)
		{
			w1_x[i][j]=w1[i][j]*Input_MinMaxNormalization[0][j];
		}
	}
	//x2 - the output of hidden layer 1 nodes
	std::vector<std::vector<double> > x2(1);
	double radbas_sum=0.0;
	for (int i = 0; i < w1.size(); i++)
	{
		double wieghtedXAndBias_sum=0.0;
		for(int j=0;j<w1[0].size();j++)
		{
			wieghtedXAndBias_sum+=w1_x[i][j];
		}
		wieghtedXAndBias_sum+=b1[i][0];
		double radbas=pow(e,-pow(wieghtedXAndBias_sum,2));
		radbas_sum+=radbas;
		x2[0].push_back(radbas);
	}
	//radbas normalization
	for (int i = 0; i < w1.size(); i++)
	{
		x2[0][i]=x2[0][i]/radbas_sum;
	}
	//w2*x2
	std::vector<std::vector<double> > w2_x(w2.size(),std::vector<double>(w2[0].size()));
	for(int i=0;i<w2.size();i++)
	{
		for(int j=0;j<w2[0].size();j++)
		{
			w2_x[i][j]=w2[i][j]*x2[0][j];
		}
	}
	//x3 - the output of hidden layer 2 nodes
	std::vector<std::vector<double> > x3(1);
	radbas_sum=0.0;
	for (int i = 0; i < w2.size(); i++)
	{
		double wieghtedXAndBias_sum=0.0;
		for(int j=0;j<w2[0].size();j++)
		{
			wieghtedXAndBias_sum+=w2_x[i][j];
		}
		wieghtedXAndBias_sum+=b2[i][0];
		double radbas=pow(e,-pow(wieghtedXAndBias_sum,2));
		radbas_sum+=radbas;
		x3[0].push_back(radbas);
	}
	//radbas normalization
	for (int i = 0; i < w2.size(); i++)
	{
		x3[0][i]=x3[0][i]/radbas_sum;
	}
	//w3*x3
	std::vector<std::vector<double> > w3_x(w3.size(),std::vector<double>(w3[0].size()));
	for(int i=0;i<w3.size();i++)
	{
		for(int j=0;j<w3[0].size();j++)
		{
			w3_x[i][j]=w3[i][j]*x3[0][j];
		}
	}
	//x4 - the output 
	std::vector<std::vector<double> > x4(1);
	for (int i = 0; i < w3.size(); i++)
	{
		double wieghtedXAndBias_sum=0.0;
		for(int j=0;j<w3[0].size();j++)
		{
			wieghtedXAndBias_sum+=w3_x[i][j];
		}
		wieghtedXAndBias_sum+=b3[i][0];
		x4[0].push_back(wieghtedXAndBias_sum);
	}
	//normalize output
	double normalizedOutput=x4[0][0];
	//reverse min max normalization
	double outputNormlized=(((normalizedOutput+1)/2)*outputMinAndRange[1][0])+outputMinAndRange[0][0];
	//reverse ave dev normalizaton
	double	output=(outputNormlized*AveDevNormalizationCoeficients[1][0])+AveDevNormalizationCoeficients[0][0];
	return output;
}
vector<vector<double>> fixturemodule::computeDefelectionArm(vector<int> indexInOutPutPath,	vector<arrayLocation>*selectionDomain,int forceIndex)
{
	string errorLevel;
	try{
		vector<double> deflectionArm_x(3,0); // from clamp 1, clamp 2, cutting force
		vector<double> deflectionArm_y(3,0); // from clamp 1, clamp 2, cutting force
		// along x
		vector<double> nearestLocatorPoint_x(3,indexInOutPutPath[0]); // initilize with locator at index 0
		errorLevel = "17.7.1";
		//Get the nearest locator point for clamps and machining force
		//along x
		//clamp1
		double loc_dist_x=abs((*selectionDomain)[indexInOutPutPath[0]].a-(*selectionDomain)[indexInOutPutPath[3]].a);
		errorLevel = "17.7.1a";
		if(loc_dist_x>abs((*selectionDomain)[indexInOutPutPath[1]].a-(*selectionDomain)[indexInOutPutPath[3]].a))
		{
			errorLevel = "17.7.b";
			loc_dist_x=abs((*selectionDomain)[indexInOutPutPath[1]].a-(*selectionDomain)[indexInOutPutPath[3]].a);
			errorLevel = "17.7.1.c";
			nearestLocatorPoint_x[0]=indexInOutPutPath[1];
		} if(loc_dist_x>abs((*selectionDomain)[indexInOutPutPath[2]].a-(*selectionDomain)[indexInOutPutPath[3]].a))
		{
			errorLevel = "17.7.1.d";
			nearestLocatorPoint_x[0]=indexInOutPutPath[2];
		}
		errorLevel = "17.7.1.e";
		//clamp2
		loc_dist_x=abs((*selectionDomain)[indexInOutPutPath[0]].a-(*selectionDomain)[indexInOutPutPath[4]].a);
		if(loc_dist_x>abs((*selectionDomain)[indexInOutPutPath[1]].a-(*selectionDomain)[indexInOutPutPath[4]].a))
		{
			loc_dist_x=abs((*selectionDomain)[indexInOutPutPath[1]].a-(*selectionDomain)[indexInOutPutPath[4]].a);
			nearestLocatorPoint_x[1]=indexInOutPutPath[1];
		} if(loc_dist_x>abs((*selectionDomain)[indexInOutPutPath[2]].a-(*selectionDomain)[indexInOutPutPath[4]].a))
		{
			nearestLocatorPoint_x[1]=indexInOutPutPath[2];
		}
		errorLevel = "17.7.1.f";
		//Cutting force
		loc_dist_x=abs((*selectionDomain)[indexInOutPutPath[0]].a*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X);
		errorLevel = "17.7.1.f.1";
		if(loc_dist_x>abs((*selectionDomain)[indexInOutPutPath[1]].a*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X))
		{
			errorLevel = "17.7.1.f.2";
			loc_dist_x=abs((*selectionDomain)[indexInOutPutPath[1]].a*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X);
			errorLevel = "17.7.1.f.3";
			nearestLocatorPoint_x[2]=indexInOutPutPath[1];
		} if(loc_dist_x>abs((*selectionDomain)[indexInOutPutPath[2]].a*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X))
		{errorLevel = "17.7.1.f.4";
		nearestLocatorPoint_x[2]=indexInOutPutPath[2];
		}
		errorLevel = "17.7.2";
		//Compute deflectin arms along x
		//clamp1 
		loc_dist_x=(double)(*selectionDomain)[nearestLocatorPoint_x[0]].a-(*selectionDomain)[indexInOutPutPath[3]].a;  	errorLevel = "17.7.2.a";
		if(loc_dist_x!=0)
		{	
			errorLevel = "17.7.2.a0";
			double alphaSteps=abs(loc_dist_x);   
			int stepChange=loc_dist_x>0?1:-1;
			errorLevel = "17.7.2.a1";
			for(int i=0;i<alphaSteps;i++)
			{
				errorLevel = "17.7.2.a2";
				deflectionArm_x[0]+=(1+3*i+3*i*i)/areaMomentOfInertia_x[(*selectionDomain)[indexInOutPutPath[3]].a+i*stepChange];
			}errorLevel = "17.7.2.a4";
		}
		//clamp2
		loc_dist_x=(double)(*selectionDomain)[nearestLocatorPoint_x[1]].a-(*selectionDomain)[indexInOutPutPath[4]].a;
		if(loc_dist_x!=0)
		{
			double alphaSteps=abs(loc_dist_x);
			int stepChange=loc_dist_x>0?1:-1;
			for(int i=0;i<alphaSteps;i++)
			{
				deflectionArm_x[1]+=(1+3*i+3*i*i)/areaMomentOfInertia_x[(*selectionDomain)[indexInOutPutPath[4]].a+i*stepChange];
			}
		}
		errorLevel = "17.7.2.b";
		//cutting force
		loc_dist_x=(double)(*selectionDomain)[nearestLocatorPoint_x[2]].a*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X;
		int alphaSteps=abs(loc_dist_x)/alpha;
		if(alphaSteps!=0)
		{
			int stepChange=loc_dist_x>0?1:-1;
			for(int i=0;i<alphaSteps;i++)
			{
				//check if outside the range
				if(areaMomentOfInertia_x.size()<=(int)real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X/alpha+i*stepChange)
				{
					UpdateResultSummary("Error in estimating defelction arm for force");
				}
				deflectionArm_x[2]+=(1+3*i+3*i*i)/areaMomentOfInertia_x[(int)real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].X/alpha+i*stepChange];
			}
		}
		//Get the nearest locator point for clamps and machining force
		//along y
		vector<double> nearestLocatorPoint_y(3,indexInOutPutPath[0]); // initilize with locator at index 0
		errorLevel = "17.7.1";
		//clamp1
		double loc_dist_y=abs((*selectionDomain)[indexInOutPutPath[0]].b-(*selectionDomain)[indexInOutPutPath[3]].b); 
		errorLevel = "17.7.1a";
		if(loc_dist_y>abs((*selectionDomain)[indexInOutPutPath[1]].b-(*selectionDomain)[indexInOutPutPath[3]].b))
		{
			errorLevel = "17.7.b";
			loc_dist_y=abs((*selectionDomain)[indexInOutPutPath[1]].b-(*selectionDomain)[indexInOutPutPath[3]].b);	
			errorLevel = "17.7.1.c";
			nearestLocatorPoint_y[0]=indexInOutPutPath[1];
		} if(loc_dist_y>abs((*selectionDomain)[indexInOutPutPath[2]].b-(*selectionDomain)[indexInOutPutPath[3]].b))
		{
			errorLevel = "17.7.1.d";
			nearestLocatorPoint_y[0]=indexInOutPutPath[2];
		}
		errorLevel = "17.7.1.e";
		//clamp2
		loc_dist_y=abs((*selectionDomain)[indexInOutPutPath[0]].b-(*selectionDomain)[indexInOutPutPath[4]].b);
		if(loc_dist_y>abs((*selectionDomain)[indexInOutPutPath[1]].b-(*selectionDomain)[indexInOutPutPath[4]].b))
		{
			loc_dist_y=abs((*selectionDomain)[indexInOutPutPath[1]].b-(*selectionDomain)[indexInOutPutPath[4]].b);
			nearestLocatorPoint_y[1]=indexInOutPutPath[1];
		} if(loc_dist_y>abs((*selectionDomain)[indexInOutPutPath[2]].b-(*selectionDomain)[indexInOutPutPath[4]].b))
		{
			nearestLocatorPoint_y[1]=indexInOutPutPath[2];
		}
		errorLevel = "17.7.1.f";
		//Cutting force
		loc_dist_y=abs((*selectionDomain)[indexInOutPutPath[0]].b*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y);  
		errorLevel = "17.7.1.f.1";
		if(loc_dist_y>abs((*selectionDomain)[indexInOutPutPath[1]].b*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y))
		{ 
			errorLevel = "17.7.1.f.2"; 
			loc_dist_y=abs((*selectionDomain)[indexInOutPutPath[1]].b*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y);
			errorLevel = "17.7.1.f.3";
			nearestLocatorPoint_y[2]=indexInOutPutPath[1];
		} if(loc_dist_y>abs((*selectionDomain)[indexInOutPutPath[2]].b*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y))
		{errorLevel = "17.7.1.f.4";
		nearestLocatorPoint_y[2]=indexInOutPutPath[2];
		}
		errorLevel = "17.7.2";
		//Compute deflectin arms along y
		loc_dist_y=(double)(*selectionDomain)[nearestLocatorPoint_y[0]].b-(*selectionDomain)[indexInOutPutPath[3]].b;
		if(loc_dist_y!=0)
		{
			double alphaSteps=abs(loc_dist_y);   
			int stepChange=loc_dist_y>0?1:-1;
			for(int i=0;i<alphaSteps;i++)
			{
				deflectionArm_y[0]+=(1+3*i+3*i*i)/areaMomentOfInertia_y[(*selectionDomain)[indexInOutPutPath[3]].b+i*stepChange];
			}
		}
		//clamp2
		loc_dist_y=(double)(*selectionDomain)[nearestLocatorPoint_y[1]].b-(*selectionDomain)[indexInOutPutPath[4]].b;
		if(loc_dist_y!=0)
		{
			double alphaSteps=abs(loc_dist_y);
			int stepChange=loc_dist_y>0?1:-1;
			for(int i=0;i<alphaSteps;i++)
			{
				deflectionArm_y[1]+=(1+3*i+3*i*i)/areaMomentOfInertia_y[(*selectionDomain)[indexInOutPutPath[4]].b+i*stepChange];
			}
		}
		//cutting force
		loc_dist_y=(double)(*selectionDomain)[nearestLocatorPoint_y[2]].b*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y;
		alphaSteps=abs(loc_dist_y)/alpha;
		if(alphaSteps!=0)
		{
			int stepChange=loc_dist_y>0?1:-1;
			for(int i=0;i<alphaSteps;i++)
			{
				//check if outside the range
				if(areaMomentOfInertia_y.size()<=(int)real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y/alpha+i*stepChange)
				{
					UpdateResultSummary("Error in estimating defelction arm for force");
				}
				deflectionArm_y[2]+=(1+3*i+3*i*i)/areaMomentOfInertia_y[(int)real_forcesApplicationPoints_fromZeroReference_Transforemed[forceIndex].Y/alpha+i*stepChange];
			}
		}
		vector<vector<double>> deflectionArm(2);
		deflectionArm[0]=deflectionArm_x;
		deflectionArm[1]=deflectionArm_y;
		return deflectionArm;
	}catch(exception ex)
	{
		messageInfo(errorLevel);
	}
}
void fixturemodule::distanceScore(int**integerOfPointsCloud,vector<	vector<arrayLocation>>* outputPath)
{
	pointsInSmallregion_all.resize((*outputPath)[0].size(),0),pointsInBigRegion_all.resize((*outputPath)[0].size(),0);
	for(int p=0;p<(*outputPath)[0].size();p++)
	{
		if(p%10==0)
		{
			double	Progress=((double)((p)/(double)((*outputPath)[0].size()))*100);
			UpdateProgressBar("Compute distanceScore progress "+to_string((int)Progress)+"%");
		}
		for (int i = 0; i < nodesNumber_H; i++) {
			for (int j = 0; j < nodesNumber_V; j++) {
				if(		integerOfPointsCloud[i][j]!=2)
				{
					double x_dif=alpha*((*outputPath)[0][p].a-i);
					double y_dif=alpha*((*outputPath)[0][p].b-j);
					double distance=sqrt(x_dif*x_dif+y_dif*y_dif);
					double intensityFactor=0.0; // weight at distance/ number of points at a distance
					if(distance<(boundaryDiagonal/2))
					{
						intensityFactor=(1-distance/(boundaryDiagonal/2));
						pointsInBigRegion_all[p]+=intensityFactor*1.0; 
					}
					if(distance<((boundaryDiagonal/2)*frictionOfSmallRegionDiamter))
					{
						intensityFactor=(1-distance/((boundaryDiagonal/2)*frictionOfSmallRegionDiamter));
						pointsInSmallregion_all[p]+=intensityFactor*1.0;
					}
				}
			}
		}
	}
	//compute distance score for forces point of action
	pointsInSmallregion_F.resize(real_forcesApplicationPoints_fromZeroReference_Transforemed.size(),0),pointsInBigRegion_F.resize(real_forcesApplicationPoints_fromZeroReference_Transforemed.size(),0);
	for(int p=0;p<real_forcesApplicationPoints_fromZeroReference_Transforemed.size();p++)
	{
		for (int i = 0; i < nodesNumber_H; i++) {
			for (int j = 0; j < nodesNumber_V; j++) {
				if(		integerOfPointsCloud[i][j]!=2)
				{
					double x_dif=(real_forcesApplicationPoints_fromZeroReference_Transforemed[p].X-i*alpha);
					double y_dif=(real_forcesApplicationPoints_fromZeroReference_Transforemed[p].Y-j*alpha);
					double distance=sqrt(x_dif*x_dif+y_dif*y_dif);
					double intensityFactor=0.0;
					//check if force is outside the geometry
					if(distance<(boundaryDiagonal/2))
					{
						intensityFactor=(1-distance/(boundaryDiagonal/2));
						pointsInBigRegion_F[p]+=intensityFactor*1;
					}
					if(distance<((boundaryDiagonal/2)*frictionOfSmallRegionDiamter))
					{
						intensityFactor=(1-distance/((boundaryDiagonal/2)*frictionOfSmallRegionDiamter));
						pointsInSmallregion_F[p]+=intensityFactor*1; 
					}
				}
			}
		}
		pointsInBigRegion_F[p]==0?pointsInBigRegion_F[p]=ebsilon:pointsInBigRegion_F[p];
		pointsInSmallregion_F[p]==0?pointsInSmallregion_F[p]=ebsilon:pointsInSmallregion_F[p];
	}
}
void fixturemodule::cavitiesScore(vector<arrayLocation> extremePoints,vector <double>  weights_x,vector <double>  weights_y,vector<	vector<arrayLocation>>* outputPath)
{
	//Cavity minimum distance and cumulative score
	minimumDistanceToCavity_all.resize((*outputPath)[0].size(),boundaryDiagonal);
	cumulativeDistanceToCavity_all.resize((*outputPath)[0].size(),0);
	for(int p=0;p<(*outputPath)[0].size();p++)
	{
		if(p%10==0)
		{
			double	Progress=((double)((p)/(double)((*outputPath)[0].size()))*100);
			UpdateProgressBar("Compute cavitiesScore progress "+to_string((int)Progress)+"%");
		}
		for(int e=numberOfextremePointsOfOuterPath;e<extremePoints.size();e++)
		{
			double x_dif=alpha*((*outputPath)[0][p].a-extremePoints[e].a);
			double y_dif=alpha*((*outputPath)[0][p].b-extremePoints[e].b);
			double distance=sqrt(x_dif*x_dif+y_dif*y_dif);
			if(minimumDistanceToCavity_all[p]==0.0)
			{
				minimumDistanceToCavity_all[p]=distance;
			}else
			{
				minimumDistanceToCavity_all[p]=min(distance,minimumDistanceToCavity_all[p]);
			}
			cumulativeDistanceToCavity_all[p]+=sqrt(weights_x[e]*weights_x[e]* x_dif*x_dif+weights_y[e]*weights_y[e]*y_dif*y_dif);
		}
		cumulativeDistanceToCavity_all[p]==0?cumulativeDistanceToCavity_all[p]=boundaryDiagonal:cumulativeDistanceToCavity_all[p];
	}
	//Forces  to Cavity minimum distance and cumulative score
	minimumDistanceToCavity_F.resize(real_forcesApplicationPoints_fromZeroReference_Transforemed.size(),boundaryDiagonal);
	cumulativeDistanceToCavity_F.resize(real_forcesApplicationPoints_fromZeroReference_Transforemed.size(),0);
	for(int p=0;p<real_forcesApplicationPoints_fromZeroReference_Transforemed.size();p++)
	{
		for(int e=numberOfextremePointsOfOuterPath;e<extremePoints.size();e++)
		{
			double x_dif=(real_forcesApplicationPoints_fromZeroReference_Transforemed[p].X-extremePoints[e].a*alpha);
			double y_dif=(real_forcesApplicationPoints_fromZeroReference_Transforemed[p].Y-extremePoints[e].b*alpha);
			double distance=sqrt(x_dif*x_dif+y_dif*y_dif);
			if( e==numberOfextremePointsOfOuterPath)
			{
				minimumDistanceToCavity_F[p]=distance;
			}else
			{
				minimumDistanceToCavity_F[p]=min(distance,minimumDistanceToCavity_F[p]);
			}
			cumulativeDistanceToCavity_F[p]+=sqrt(weights_x[e]*weights_x[e]* x_dif*x_dif+weights_y[e]*weights_y[e]*y_dif*y_dif);
		}
		cumulativeDistanceToCavity_F[p]==0?cumulativeDistanceToCavity_F[p]=boundaryDiagonal:cumulativeDistanceToCavity_F[p];
	}
}
void fixturemodule::conicalBeamAndShearLine(int**integerOfPointsCloud,vector<	vector<arrayLocation>> *outputPath)
{
	//Points within conical beam
	pointsInForawrdConicalBeam_all.resize((*outputPath)[0].size(),0);
	lengthOfShearLine_all.resize((*outputPath)[0].size(),0);
	for(int p=0;p<(*outputPath)[0].size();p++)
	{
		if(p%10==0)
		{
			double	Progress=((double)((p)/(double)((*outputPath)[0].size()))*100);
			UpdateProgressBar("Compute points in forward beam progress "+to_string((int)Progress)+"%");
		}
		double distanceToNearestCavityWithinBeam=0.0;
		double distanceToFurthestPointWithinBeam=0.0;
		double firstOnSurfaceCollinerPoint_trigger=false;
		for (int i = 0; i < nodesNumber_H; i++) {
			for (int j = 0; j < nodesNumber_V; j++) {
				vector<double> firstVector(2),secondVector(2);
				firstVector[0]=pathPointsUnitNorm_x_all[0][p];
				firstVector[1]=pathPointsUnitNorm_y_all[0][p];
				secondVector[0]=alpha*( (*outputPath)[0][p].a-i);
				secondVector[1]=alpha*( (*outputPath)[0][p].b-j);
				double dot= firstVector[0]*secondVector[0]+firstVector[1]*secondVector[1]; // x1*x2 + y1*y2
				double det= firstVector[0]*secondVector[1]-firstVector[1]*secondVector[0]; // x1*y2 - y1*x2
				double angleInRadian=atan2(det,dot);
				double lengthOfSecondVector=sqrt(  secondVector[0]*secondVector[0]+  secondVector[1]*secondVector[1]   );
				double intensityFactor=(1-lengthOfSecondVector/boundaryDiagonal);
				if(		integerOfPointsCloud[i][j]!=2)
				{
					if(abs(angleInRadian)<coneAngleInRadian)
					{
						pointsInForawrdConicalBeam_all[p]+=1*intensityFactor;
						if(abs(tan(angleInRadian))*lengthOfSecondVector<alpha)
						{
							distanceToFurthestPointWithinBeam=max(distanceToFurthestPointWithinBeam,lengthOfSecondVector);
							firstOnSurfaceCollinerPoint_trigger=true;
						}
					}
				}else
				{
					if(abs(angleInRadian)<coneAngleInRadian)
					{
						if(abs(tan(angleInRadian))*lengthOfSecondVector<alpha)
						{
							if(lengthOfSecondVector!=0 && firstOnSurfaceCollinerPoint_trigger)
							{
								if(distanceToNearestCavityWithinBeam==0)
								{
									distanceToNearestCavityWithinBeam=lengthOfSecondVector;
								}else
								{
									distanceToNearestCavityWithinBeam=min(distanceToNearestCavityWithinBeam,lengthOfSecondVector);
								}	
							}	
						}
					}
				}
			}
		}
		if(distanceToNearestCavityWithinBeam!=0)
		{
			lengthOfShearLine_all[p]=distanceToNearestCavityWithinBeam;
		}else
		{
			lengthOfShearLine_all[p]=distanceToFurthestPointWithinBeam;
		}
	}
	//(for forces) Points within conical beam
	pointsInForawrdConicalBeam_F.resize(real_forcesApplicationPoints_fromZeroReference_Transforemed.size(),0);
	lengthOfShearLine_F.resize(real_forcesApplicationPoints_fromZeroReference_Transforemed.size(),0);
	for(int p=0;p<real_forcesApplicationPoints_fromZeroReference_Transforemed.size();p++)
	{
		double distanceToNearestCavityWithinBeam=0.0;
		double distanceToFurthestPointWithinBeam=0.0;
		double firstOnSurfaceCollinerPoint_trigger=false;
		for (int i = 0; i < nodesNumber_H; i++) {
			for (int j = 0; j < nodesNumber_V; j++) {
				vector<double> firstVector(2),secondVector(2);
				firstVector[0]=appliedForcesMoments_atZeroRef[p][0][0];
				firstVector[1]=appliedForcesMoments_atZeroRef[p][1][0];
				secondVector[0]= i*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[p].X;
				secondVector[1]= j*alpha-real_forcesApplicationPoints_fromZeroReference_Transforemed[p].Y;
				double dot= firstVector[0]*secondVector[0]+firstVector[1]*secondVector[1]; // x1*x2 + y1*y2
				double det= firstVector[0]*secondVector[1]-firstVector[1]*secondVector[0]; // x1*y2 - y1*x2
				double angleInRadian=atan2(det,dot);
				double lengthOfSecondVector=sqrt(  secondVector[0]*secondVector[0]+  secondVector[1]*secondVector[1]   );
				double intensityFactor=(1-lengthOfSecondVector/boundaryDiagonal);
				if(		integerOfPointsCloud[i][j]!=2)
				{
					if(abs(angleInRadian)<coneAngleInRadian)
					{
						pointsInForawrdConicalBeam_F[p]+=1*intensityFactor;
						if(abs(tan(angleInRadian))*lengthOfSecondVector<alpha)
						{
							distanceToFurthestPointWithinBeam=max(distanceToFurthestPointWithinBeam,lengthOfSecondVector);
							firstOnSurfaceCollinerPoint_trigger=true;
						}
					}
				}else
				{
					if(abs(angleInRadian)<coneAngleInRadian)
					{
						if(abs(tan(angleInRadian))*lengthOfSecondVector<alpha)
						{
							if(lengthOfSecondVector!=0 && firstOnSurfaceCollinerPoint_trigger)
							{
								if(distanceToNearestCavityWithinBeam==0)
								{
									distanceToNearestCavityWithinBeam=lengthOfSecondVector;
								}else
								{
									distanceToNearestCavityWithinBeam=min(distanceToNearestCavityWithinBeam,lengthOfSecondVector);
								}	
							}	
						}
					}
				}
			}
		}
		if(distanceToNearestCavityWithinBeam!=0)
		{
			lengthOfShearLine_F[p]=distanceToNearestCavityWithinBeam;
		}else
		{
			lengthOfShearLine_F[p]=distanceToFurthestPointWithinBeam;
		}
	}
}
vector<std::vector<double>> fixturemodule::annInputs(vector<int> generatedPointsLocation,	int pointsInOuterEdge)
{
	string errorLevel;
	try{
		vector<std::vector<double>>  inputx;
		for (int p = 0; p < limda1.size(); p++)
		{
			vector<double> input;
			input.push_back(boundaryDiagonal/alpha); //1 
			double cuttingForceToClampingForceRatio=sqrt(appliedForcesMoments_atZeroRef[p][0][0]*appliedForcesMoments_atZeroRef[p][0][0]+
				appliedForcesMoments_atZeroRef[p][1][0]*appliedForcesMoments_atZeroRef[p][1][0])/clampingForce;
			input.push_back(cuttingForceToClampingForceRatio); //2
			input.push_back(limda1[p]/clampingForce);  //3
			input.push_back(limda2[p]/clampingForce);  //4
			input.push_back(limda3[p]/clampingForce);  //5
			input.push_back(cumulativeCrossDistanceOfFixels/boundaryDiagonal); //6
			for (int i = 0; i < 5; i++)
			{
				input.push_back(pointsInSmallregion_subSet[generatedPointsLocation[i]]); //7,8,9,10,11
			}
			input.push_back(pointsInSmallregion_F[p]); //12
			for (int i = 0; i < 5; i++)
			{
				input.push_back(pointsInBigRegion_subSet[generatedPointsLocation[i]]); //13,14,15,16,17
			}
			input.push_back(pointsInBigRegion_F[p]); //18
			for (int i = 0; i < 5; i++)
			{
				input.push_back(pointsInForawrdConicalBeam_subSet[generatedPointsLocation[i]]); //19,20,21,22,23
			}
			input.push_back(pointsInForawrdConicalBeam_F[p]); //24
			for (int i = 0; i < 5; i++)
			{
				input.push_back(lengthOfShearLine_subSet[generatedPointsLocation[i]]/boundaryDiagonal); //25,26,27,28,29
			}
			input.push_back(lengthOfShearLine_F[p]/boundaryDiagonal); //30
			vector<double> freedomScore(3,0);
			//Clamp 1 freedom - Cumulative distance from the three locators
			for (int i = 0; i < 3; i++)
			{
				double dist_x=alpha*(selectionDomain[generatedPointsLocation[3]].a-selectionDomain[generatedPointsLocation[i]].a);
				double dist_y=alpha*(selectionDomain[generatedPointsLocation[3]].b-selectionDomain[generatedPointsLocation[i]].b);
				freedomScore[0]+=sqrt(dist_x*dist_x+dist_y*dist_y);
			}
			//Clamp 2 freedom - Cumulative distance from the three locators
			for (int i = 0; i < 3; i++)
			{
				double dist_x=alpha*(selectionDomain[generatedPointsLocation[4]].a-selectionDomain[generatedPointsLocation[i]].a);
				double dist_y=alpha*(selectionDomain[generatedPointsLocation[4]].b-selectionDomain[generatedPointsLocation[i]].b);
				freedomScore[1]+=sqrt(dist_x*dist_x+dist_y*dist_y);
			}
			//Cutting force freedom - Cumulative distance from the three locators
			for (int i = 0; i < 3; i++)
			{
				double dist_x=real_forcesApplicationPoints_fromZeroReference_Transforemed[p].X-alpha*selectionDomain[generatedPointsLocation[i]].a;
				double dist_y=real_forcesApplicationPoints_fromZeroReference_Transforemed[p].Y-alpha*selectionDomain[generatedPointsLocation[i]].b;
				freedomScore[2]+=sqrt(dist_x*dist_x+dist_y*dist_y);
			}
			for (int i = 0; i < 3; i++)
			{
				input.push_back(freedomScore[i]/boundaryDiagonal); //31,32,33
			}
			for (int i = 0; i < 3; i++)
			{
				double deflectionResultant=sqrt((deflectionScore_x[p][i]*deflectionScore_x[p][i])+(deflectionScore_y[p][i]*deflectionScore_y[p][i]));
				input.push_back(deflectionResultant/boundaryDiagonal);	//34.35,36
			}
			saveStringToFile(&cal.matrixToString(input),"annInput");
			inputx.push_back(input);
		}
		return inputx;
	}catch(exception ex)
	{ 
		messageInfo(errorLevel);
	}
}

void fixturemodule::fixelsDistanceAttributes(vector<int> generatedPointsLocation,	vector<arrayLocation>* selectionDomain)
{
	// cross attributes
	vector<double> fixelsCrossDistances;
	for(int i=0;i<5;i++)
	{
		for(int j=i+1;j<5;j++)
		{
			double x_dist=alpha*(double)( (*selectionDomain)[generatedPointsLocation[i]].a-(*selectionDomain)[generatedPointsLocation[j]].a);
			double y_dist=alpha*(double)( (*selectionDomain)[generatedPointsLocation[i]].b-(*selectionDomain)[generatedPointsLocation[j]].b);
			double dist=sqrt(x_dist*x_dist+y_dist*y_dist);
			fixelsCrossDistances.push_back(dist);
		}
	}
	cumulativeCrossDistanceOfFixels=accumulate(std::begin(fixelsCrossDistances), std::end(fixelsCrossDistances),0 );
	minCrossDistanceOfFixels= *min_element(std::begin(fixelsCrossDistances), std::end(fixelsCrossDistances) );
	maxCrossDistanceOfFixels=*max_element(std::begin(fixelsCrossDistances), std::end(fixelsCrossDistances) );
	distanceDeviationcross=cal.standardDeviation(fixelsCrossDistances);
	//In path attributes
	vector< int > tempOrder=generatedPointsLocation;
	sort(tempOrder.begin(),tempOrder.end());
	vector<double> fixelsAlongPathDistances;
	for(int i=0;i<5;i++)
	{
		double x_dist=alpha*(double)( (*selectionDomain)[tempOrder[i]].a-(*selectionDomain)[tempOrder[(i+1)%5]].a);
		double y_dist=alpha*(double)( (*selectionDomain)[tempOrder[i]].b-(*selectionDomain)[tempOrder[(i+1)%5]].b);
		double dist=sqrt(x_dist*x_dist+y_dist*y_dist);
		fixelsAlongPathDistances.push_back(dist);
	}
	cumulativeInPathDistanceOfFixels=accumulate(std::begin(fixelsAlongPathDistances), std::end(fixelsAlongPathDistances),0 );
	minInPathDistanceOfFixels=*min_element(std::begin(fixelsAlongPathDistances), std::end(fixelsAlongPathDistances) );
	maxInPathDistanceOfFixels=*max_element(std::begin(fixelsAlongPathDistances), std::end(fixelsAlongPathDistances) );
	distanceDeviationAlongPath=cal.standardDeviation(fixelsAlongPathDistances);
}
void fixturemodule::deflictionScore(vector<int> generatedPointsLocation)
{
	deflectionScore.resize(appliedForcesMoments_atZeroRef.size(),vector<double>(3));
	deflectionScore_x.resize(appliedForcesMoments_atZeroRef.size(),vector<double>(3));
	deflectionScore_y.resize(appliedForcesMoments_atZeroRef.size(),vector<double>(3));
	for (int i = 0; i < appliedForcesMoments_atZeroRef.size(); i++)
	{
		//defliction score preperation
		double forceMagnitude=sqrt(appliedForcesMoments_atZeroRef[i][0][0]*appliedForcesMoments_atZeroRef[i][0][0]+
			appliedForcesMoments_atZeroRef[i][1][0]*appliedForcesMoments_atZeroRef[i][1][0]);
		//	vector<double> deflectionScore_x(3,0); 
		deflectionScore_x[i][0]= pathPointsUnitNorm_y_subSet[generatedPointsLocation[3]]*alpha*alpha*alpha*deflectionArms[i][0][0];
		deflectionScore_x[i][1]=pathPointsUnitNorm_y_subSet[generatedPointsLocation[4]]*alpha*alpha*alpha*deflectionArms[i][0][1];
		if(forceMagnitude!=0)
		{
			deflectionScore_x[i][2]=(appliedForcesMoments_atZeroRef[i][1][0]/forceMagnitude)*alpha*alpha*alpha*deflectionArms[i][0][2];
		}
		//	vector<double> deflectionScore_y(3,0);
		deflectionScore_y[i][0]= pathPointsUnitNorm_x_subSet[generatedPointsLocation[3]]*alpha*alpha*alpha* deflectionArms[i][1][0];
		deflectionScore_y[i][1]= pathPointsUnitNorm_x_subSet[generatedPointsLocation[4]]*alpha*alpha*alpha* deflectionArms[i][1][1];
		if(forceMagnitude!=0)
		{
			deflectionScore_y[i][2]=(appliedForcesMoments_atZeroRef[i][0][0]/forceMagnitude)*alpha*alpha*alpha* deflectionArms[i][1][2];
		}
		/*deflectionScore[i][0]=sqrt(deflectionScore_x[0]*deflectionScore_x[0]+deflectionScore_y[0]*deflectionScore_y[0]);
		deflectionScore[i][1]=sqrt(deflectionScore_x[1]*deflectionScore_x[1]+deflectionScore_y[1]*deflectionScore_y[1]);
		deflectionScore[i][2]=sqrt(deflectionScore_x[2]*deflectionScore_x[2]+deflectionScore_y[2]*deflectionScore_y[2]);*/
	}
}
void fixturemodule::saveStringToFile(string * text,string fileName)
{
	std::ofstream test("\\DeltaFix tool files\\output\\"+fileName+".txt");
	test <<*text;
	test.close();
}
void fixturemodule::showFinalResult(tag_t** pointsCloud,	vector<int> optimumPointsLocation)
{
	string errorLevel;
	vector<	Point3d> locations_transformed(5);
	vector<	Point3d> locations_orginal(5);
	UpdateResultSummary("Indexes of optimum points: "+cal.matrixToOneLine(optimumPointsLocation));
	for (int i = 0; i < locations_transformed.size(); i++)
	{
		locations_transformed[i].X=zeroLocation_transformed.X+(double)selectionDomain[optimumPointsLocation[i]].a*alpha;
		locations_transformed[i].Y=zeroLocation_transformed.Y+(double)selectionDomain[optimumPointsLocation[i]].b*alpha;
		locations_transformed[i].Z=transformed_z_distance;
		locations_orginal[i]=ReversedPoint(locations_transformed[i]);
	}
	UpdateResultSummary_seprator();
	result.push_back("Optimum locations of locators and clamps");
	result.push_back(("  L 1: ("+to_string( locations_transformed[0].X)+","+to_string( locations_transformed[0].Y)+","+to_string( locations_transformed[0].Z)+")").c_str());
	result.push_back(("  L 2: ("+to_string( locations_transformed[1].X)+","+to_string( locations_transformed[1].Y)+","+to_string( locations_transformed[1].Z)+")").c_str());
	result.push_back(("  L 3: ("+to_string( locations_transformed[2].X)+","+to_string( locations_transformed[2].Y)+","+to_string( locations_transformed[2].Z)+")").c_str());
	result.push_back((	" C 1: ("+to_string( locations_transformed[3].X)+","+to_string( locations_transformed[3].Y)+","+to_string( locations_transformed[3].Z)+")").c_str());
	result.push_back(("  C 2: ("+to_string( locations_transformed[4].X)+","+to_string( locations_transformed[4].Y)+","+to_string( locations_transformed[4].Z)+")").c_str());
	result.push_back("Elapsed time: "+elapsedTime);
	UpdateResultSummary();
	//export result to txt file
	string summary="";
	for (int i = 0; i < result.size(); i++)
	{
		summary+=result[i].GetText();
		summary+="\n";
	}
	saveStringToFile(&summary,"Result summary");
}
vector<int>  fixturemodule::DNSAExecution(vector<arrayLocation> *extremePoints,vector <double>*curveture_subSet,vector <double>  *weights_x,vector <double> *weights_y,int pointsInOuterEdge)
{
	UpdateProgressBar("DNSA execution begin");
	string errorLevel;
	try{
		int currentHighlightedGroup[5]={0};
		vector<int> optimumPointsLocation(5,0);
		vector<	vector<double>> obj_calibration(6);
		vector<double> obj_ave(6,0);
		vector<	double> obj_dev(6,1);
		double optimizationProgress=0;
		int domainSize=selectionDomain.size();
		UpdateResultSummary("Domain size = "+to_string(domainSize));
		UpdateResultSummary("Number of reference points = "+to_string(extremePoints->size()));

		for (int i = 0; i < 6; i++)
		{

			UpdateResultSummary("k"+to_string(i+1)+" = "+cal.doubleToString(k_obj[i]));

		}

		for(int c=0;c<2;c++)
		{
			if(c==0)
			{
				UpdateProgressBar("Callibartion phase of DNSA algorithm");
				UpdateResultSummary_seprator();
				UpdateResultSummary("Callibartion phase started");
				UpdateResultSummary_seprator();
			}else
			{
				optimizationProgress=0;
				UpdateProgressBar("Optimization progress "+to_string((int)optimizationProgress)+"%");
			}
			double minimumCost=0;
			vector<double> pointMeanPointer(5,0);
			vector<vector<int>> optimumPointsSetMap(5);
			vector<double> costReciprocalOfOptimum;
			//double minCost=0,maxCost=0,aveCost=0,aveCost_couter=0;
			double minCostChange=0;
			errorLevel=3;
			vector<	int> generatedPointsLocation(5,0);
			vector<int> newPoint(5,0);
			double costChange;
			double criterian;
			bool condition;
			double layoutCost;
			localizationCost=0;
			detachmentScore=0;
			locatorsReactions=0;
			double pointsSpreed;


			if(c==1)
			{
				for (int i = 0; i < obj_ave.size(); i++)
				{
					obj_ave[i]=cal.average(obj_calibration[i]);
					obj_dev[i]=cal.standardDeviation(obj_calibration[i]);
					UpdateResultSummary("obj"+to_string(i+1)+" Ave. = "+cal.doubleToString(obj_ave[i]));
					UpdateResultSummary("obj"+to_string(i+1)+" SD. = "+cal.doubleToString(obj_dev[i]));
				}
				UpdateResultSummary("Main optimization phase started");
			}


			vector<string> meansRecord(5);
			vector<string> meansRecord_300(5);
			vector<string> iterationsRecord(5);
			vector<string> zscoreRecords(6);
			vector<string> subCostsRecords(6);
			vector<string> zscoreRecords_300(6);
			string metroplisCreterian;
			string metroplisCreterian_300;
			string metroplisCreterian2;
			string metroplisCreterian2_300;
			string costRecord;
			string costRecord_300;
			string acceptedCost;
			vector<string> OptimumZscore(6);
			vector<string> OptimumCosts(6);
			string sizeOfGuideGroup;
			int maximumIterations=1000000;
			int iterationsCounter=0;



			sizeOfRandomPhase=(iterations*epochs)*precentageOfRandomSize;
			for(int d=0;d<iterations;d++)
			{
				for(int f=0;f<epochs;f++)
				{
					if(iterationsCounter==maximumIterations)
					{
						goto terminateSA;
					}
					iterationsCounter+=1;
					errorLevel="d1";
					//accelerate at callibration phase
					if(c==0)
					{
						int actualPhaseSizetoCalibratioPhase=10;
						if((d*epochs+f)%actualPhaseSizetoCalibratioPhase!=0)
						{
							continue;
						}
					}
generateNew:
					localizationCost=0; detachmentScore=0; locatorsReactions=0;limdas_spreed=0;
					double k_cooling=log10(T_max/T_min)/((double)iterations*(1-(precentageOfRandomSize/(double)epochs))-1);
					double			k=pow(e,-(k_cooling*((double)(d-((double)iterations*precentageOfRandomSize)))));
					errorLevel="d2";
					//if(f==0){	UpdateResultSummary("k="+to_string(k));}
					bool isRandomPhase;
					if((d*epochs+f)<sizeOfRandomPhase )
					{
						for (int i=0;i<5;i++)
						{
re:

							if(i==0 )
							{
								generatedPointsLocation[i]=(int)cal.randMToN(0,(double)domainSize-0.001);

							}else
							{
								generatedPointsLocation[i]=((int)generatedPointsLocation[i-1]+(int)((double)domainSize/cal.randMToN(2,3)))%(domainSize);
							}
							isRandomPhase=true;
							/*generatedPointsLocation[i]=(int)cal.randMToN(0,(double)domainSize-0.001);
							if(c==1)
							{
							if((i==0 ||i==1)&& !(selectionDomain[generatedPointsLocation[i]].a<nodesNumber_H-3 && selectionDomain[generatedPointsLocation[i]].b==0))
							{
							goto re;
							}else if((i==2)&& !(selectionDomain[generatedPointsLocation[i]].b<nodesNumber_V-3 && selectionDomain[generatedPointsLocation[i]].a==0))
							{
							goto re;
							}else if((i==3)&& !(selectionDomain[generatedPointsLocation[i]].b<nodesNumber_V-3 && selectionDomain[generatedPointsLocation[i]].a==nodesNumber_H-1))
							{
							goto re;
							}else if((i==4)&& !(selectionDomain[generatedPointsLocation[i]].a<nodesNumber_H-3 && selectionDomain[generatedPointsLocation[i]].b==nodesNumber_V-1))
							{
							goto re;
							}
							}*/
						}
					}else
					{
						//Declined PDF phase
						for (int i=0;i<5;i++)
						{
							newPoint[i]=	cal.newPointGenerator(k,(int)pointMeanPointer[i],domainSize-1);
							errorLevel="d3";
re1:		;
							//constrin
							/*if(c==1)
							{
							if((i==0 ||i==1)&& !(selectionDomain[newPoint[i]].a<nodesNumber_H-3 && selectionDomain[newPoint[i]].b==0))
							{newPoint[i]=(int)cal.randMToN(0,domainSize-1);
							goto re1;
							}else if((i==2)&& !(selectionDomain[newPoint[i]].b<nodesNumber_V-3 && selectionDomain[newPoint[i]].a==0))
							{newPoint[i]=(int)cal.randMToN(0,domainSize-1);
							goto re1;
							}else if((i==3)&& !(selectionDomain[newPoint[i]].b<nodesNumber_V-3 && selectionDomain[newPoint[i]].a==nodesNumber_H-1))
							{newPoint[i]=(int)cal.randMToN(0,domainSize-1);
							goto re1;
							}else if((i==4)&& !(selectionDomain[newPoint[i]].a<nodesNumber_H-3 && selectionDomain[newPoint[i]].b==nodesNumber_V-1))
							{newPoint[i]=(int)cal.randMToN(0,domainSize-1);
							goto re1;
							}
							errorLevel="d4";
							}*/
						}
						for (int i=0;i<newPoint.size();i++)
						{
							generatedPointsLocation[i]=newPoint[i];
						}
						isRandomPhase=false;
					}
					errorLevel="x2";
					double neuralNetworkoutput;
					//identical values check
					bool checkDublication=false;
					for(int i=0;i<  5;i++)
					{
						for(int j=i+1;j<  5;j++)
						{
							if(generatedPointsLocation[i]==generatedPointsLocation[j])
							{
								checkDublication=true;
							}
						}
					}
					if(checkDublication)
					{
						f-=1;
						continue;
					}
					bool checkOutSideRange=false;
					for(int i=0;i<  newPoint.size();i++)
					{
						if(generatedPointsLocation[i]<0 || generatedPointsLocation[i]>=domainSize)
						{
							checkOutSideRange=true;
							UpdateResultSummary("Fatal error: Generated point outside the range of the selection domain, #iteration="+to_string(d)+" #epoch="+to_string(f));
							UpdateResultSummary("generated point index: "+to_string(generatedPointsLocation[i]));
							UpdateResultSummary("generated point for fixel: "+to_string(i));
							UpdateResultSummary("constant k: "+to_string(k));
							UpdateResultSummary("domain size: "+to_string(domainSize));
							UpdateResultSummary("points Mean index"+cal.matrixToOneLine( pointMeanPointer));
							for (int i = 0; i < optimumPointsSetMap[0].size(); i++)
							{
								string pointsIndexes="";
								for (int j = 0; j < 5; j++)
								{
									pointsIndexes+=to_string(optimumPointsSetMap[j][i])+" , ";
								}
								UpdateResultSummary("Container layout "+to_string(i)+" indexes->"+pointsIndexes);
								UpdateResultSummary("Container layout "+to_string(i)+" fittness="+cal.doubleToString( costReciprocalOfOptimum[i]));
							}
							if(isRandomPhase)
							{
								UpdateResultSummary("At random phase");
							}else
							{
								UpdateResultSummary("At declinig phase");
							}
							DeleteAllPoints(pointsCloud2_sub);
							return optimumPointsLocation;
						}
					}
					if(checkOutSideRange)
					{
						f-=1;
						continue;
					}

					vector<	double> zScore(6);
					errorLevel="x2.1 d="+to_string(d)+"  f="+to_string(f)+" c="+to_string(c);
					if(	costCalculation(&pathPointsUnitNorm_x_subSet,&pathPointsUnitNorm_y_subSet,generatedPointsLocation,extremePoints,&selectionDomain,curveture_subSet,weights_x,weights_y)==0)
					{
						f-=1;
						continue;
					}
					if(c==1 && showSimulation  && (d*epochs+f)%(70)==0)
					{
						for(int i=0;i<  newPoint.size();i++)
						{
							if(currentHighlightedGroup[i]>=0 &&currentHighlightedGroup[i]<pointsCloud2_sub.size())
							{
								highlightObj(pointsCloud2_sub[currentHighlightedGroup[i]],false);
								highlightObj(pointsCloud2_sub[generatedPointsLocation[i]],true);
								currentHighlightedGroup[i]=generatedPointsLocation[i];
								UF_initialize();
								UF_MODL_update();
								UF_DISP_refresh();
								UF_terminate();
							}
						}
					}
					if(((d*epochs)+f)%(int)((iterations*epochs)/100)==0 && c==1)
					{
						optimizationProgress=((double)((d*epochs)+f+1)/(double)(iterations*epochs))*100;
						UpdateProgressBar("Optimization progress "+to_string((int)optimizationProgress)+"% ,i="+to_string(d));
					}else if(((d*epochs)+f)%(int)((iterations*epochs)/100)==0 && c==0)
					{
						optimizationProgress=((double)((d*epochs)+f+1)/(double)(iterations*epochs))*100;
						UpdateProgressBar("Callibration progress "+to_string((int)optimizationProgress)+"%");
					}
					if(localizationCost==0.0)
					{
						f-=1;
						continue;
					}
					errorLevel="x2.1.a ";
					//contact points spreedness
					vector<double> fixlesDist;
					pointsSpreed=0;
					for(int i=0;i<   5;i++)
					{
						for(int j=i+1;j<  5;j++)
						{
							double distd=abs((double)generatedPointsLocation[i]-(double)generatedPointsLocation[j]);
							distd=min(distd,(double)domainSize-distd);
							pointsSpreed+=(((double)domainSize/distd)*((double)domainSize/distd))/100;
							//fixlesDist.push_back((double)generatedPointsLocation[i]-(double)generatedPointsLocation[j]);
						}
					}
					//pointsSpreed=	1/cal.standardDeviation(fixlesDist);
					///////////////////////////////////////////////////////////////////////////////
					//prepare other attributes for ANN
					if(k_obj[2]==0.0)
					{
						neuralNetworkoutput=1;
					}else
					{
						fixelsDistanceAttributes( generatedPointsLocation,&selectionDomain);
						errorLevel="x2.1.b ";
						//defelction arm for all machining force state
						for (int i = 0; i < real_forcesApplicationPoints_fromZeroReference_Transforemed.size(); i++)
						{
							deflectionArms.push_back(computeDefelectionArm(generatedPointsLocation,&selectionDomain,i));
						}
						errorLevel="x2.1.c   israndomphase=" +to_string( isRandomPhase)+" "+ cal.matrixToString(generatedPointsLocation)+" e"+cal.matrixToString(pointMeanPointer);
						deflictionScore( generatedPointsLocation);
						vector<	vector<double>> input;
						input=	annInputs(generatedPointsLocation,pointsInOuterEdge);
						errorLevel="x2.1.c   3";
						neuralNetworkoutput=0;
						for (int i = 0; i < input.size(); i++)
						{
							vector<		vector<double>> inputx(1);
							inputx[0]=input[i]; 
							//Avoid minius values of neural network output by adding a threshold
							neuralNetworkoutput+=exposureFriction[i]*annOutput2(&inputx)+pow(10,-6);
						}
						if(neuralNetworkoutput!=neuralNetworkoutput)
						{
							UpdateResultSummary("Fatal error: undetermined neural network output");
							DeleteAllPoints(pointsCloud2_sub);
							continue;
							//return optimumPointsLocation;
						}

					}
					////////////////////////////////////////////////////////////////////////////////////////////////////////
					errorLevel="x2.1.c   4";
					/*zScore[0]=(localizationCost)/obj_dev[0];
					zScore[1]=(detachmentScore+1)/obj_dev[1];
					zScore[2]=(neuralNetworkoutput)/obj_dev[2];
					zScore[3]=(locatorsReactions+1)/obj_dev[3];
					zScore[4]=(limdas_spreed)/obj_dev[4];
					*/
					zScore[0]=(localizationCost-obj_ave[0])/obj_dev[0];
					zScore[1]=(detachmentScore-obj_ave[1])/obj_dev[1];
					if(k_obj[2]==0.0)
					{	zScore[2]=1;}else
					{
						zScore[2]=(neuralNetworkoutput-obj_ave[2])/obj_dev[2];
					}
					zScore[3]=(locatorsReactions-obj_ave[3])/obj_dev[3];
					zScore[4]=(limdas_spreed-obj_ave[4])/obj_dev[4];
					zScore[5]= (pointsSpreed-obj_ave[5])/obj_dev[5];
					//check high zscore
					bool highZscoreFound=false;
					if(c==1)
					{
						for (int i = 0; i < zScore.size(); i++)
						{
							if(zScore[i]>highZscoreTrigger)
							{
								highZscoreFound=true;
							}
						}
					}
					if(highZscoreFound)
					{
						f+=1;
						continue;
					}
					if(c==1){
						for (int i = 0; i < 6; i++)
						{
							zscoreRecords[i]+=to_string(zScore[i])+"\n";
							if(f==0)
							{
								zscoreRecords_300[i]+=to_string(zScore[i])+"\n";
							}
						}
						subCostsRecords[0]+=cal.doubleToString(localizationCost)+"\n";
						subCostsRecords[1]+=cal.doubleToString(detachmentScore)+"\n";
						subCostsRecords[2]+=cal.doubleToString(neuralNetworkoutput)+"\n";
						subCostsRecords[3]+=cal.doubleToString(locatorsReactions)+"\n";
						subCostsRecords[4]+=cal.doubleToString(limdas_spreed)+"\n";
						subCostsRecords[5]+=cal.doubleToString(pointsSpreed)+"\n";
					}

					for (int i = 0; i < zScore.size(); i++)
					{
						if(zScore[i]!=zScore[i])
						{
							UpdateResultSummary("Fatal error: Undeterninenet zscore value at i="+to_string(i+1));
						}
						/*if(zScore[i]<0.0)
						{
						UpdateResultSummary("Fatal error: Negative objective"+to_string(i+1)+" value ="+cal.doubleToString(zScore[i]));
						}else if(zScore[i]==0.0)
						{
						UpdateResultSummary("Fatal error: Zero objective"+to_string(i+1)+" value ="+cal.doubleToString(zScore[i]));
						}*/
					}
					if(c==1)
					{		

						//layoutCost=1;
						//for (int i = 0; i < zScore.size(); i++)
						//{
						//	//layoutCost+=zScore[i]*k_obj[i];
						//	layoutCost*=pow(e,(zScore[i]*k_obj[i]));

						//}
						layoutCost=0;
						for (int i = 0; i < zScore.size(); i++)
						{
							//layoutCost+=zScore[i]*k_obj[i];
							layoutCost+=pow(pow(e,zScore[i])*k_obj[i],2);
						}
						layoutCost=sqrt(layoutCost);
						if(f==0)
						{
							costRecord_300+=cal.doubleToString(layoutCost)+"\n";
						}
						costRecord+=cal.doubleToString(layoutCost)+"\n";
					}else
					{	
						/*for (int i = 0; i < 5; i++)
						{
						UpdateResultSummary("zscore "+to_string(i)+":"+ cal.doubleToString(zScore[i]));
						}*/
						layoutCost=pow(zScore[0],k_obj[0])*pow( zScore[1]+1,k_obj[1])*pow(zScore[2],k_obj[2])*pow( zScore[3]+1,k_obj[3])*pow( zScore[4]+1,k_obj[4])*pow( zScore[5]+1,k_obj[5]);
						//layoutCost=zScore[0]*zScore[1]*zScore[2]*zScore[3]*zScore[4];

					}
					//UpdateResultSummary(to_string( d));
					//calibration data
					if(c==0)
					{
						for (int i = 0; i < zScore.size(); i++)
						{
							obj_calibration[i].push_back( zScore[i]);
						}
					}
					errorLevel="x2.1.3";
					costChange=abs(minimumCost-layoutCost);
					errorLevel="x2.1.4";
					if(costChange!=0 )
					{
						if(minCostChange==0 || minCostChange>costChange)
						{
							minCostChange=costChange;
						}
						/*if(maxCost<costChange)
						{
						maxCost=costChange;
						}
						aveCost=((aveCost*aveCost_couter )+costChange)/(aveCost_couter+1) ;
						aveCost_couter+=1;*/
						errorLevel="x2.1.5";
					}
					bool newacceptedCandidate=false;
					if( minimumCost>layoutCost ||minimumCost==0)
					{
						if(c==1){
							metroplisCreterian2+=to_string(1)+"\n";
							if(f==0)
							{
								metroplisCreterian2_300+=to_string(1)+"\n";
							}
							iterationsRecord[0]+=to_string(d)+"\n";
							for (int i = 0; i < 6; i++)
							{
								OptimumZscore[i]=cal.doubleToString( zScore[i]);
							}
							OptimumCosts[0]=cal.doubleToString(localizationCost);
							OptimumCosts[1]=cal.doubleToString(detachmentScore);
							OptimumCosts[2]=cal.doubleToString(neuralNetworkoutput);
							OptimumCosts[3]=cal.doubleToString(locatorsReactions);
							OptimumCosts[4]=cal.doubleToString(limdas_spreed);
							OptimumCosts[5]=cal.doubleToString(pointsSpreed);
						}
						errorLevel="x4.n";
						//تحديث الموقع الامثل
						minimumCost=layoutCost; 
						for(int i=0;i<  5;i++)
						{
							optimumPointsLocation[i]=generatedPointsLocation[i];
						}
						limda1_optimum.clear();limda2_optimum.clear();limda3_optimum.clear();
						for (int i = 0; i < limda1.size(); i++)
						{
							limda1_optimum.push_back(limda1[i]);
							limda2_optimum.push_back(limda2[i]);
							limda3_optimum.push_back(limda3[i]);
						}
						errorLevel="x4.n1";
						for(int i=0;i<  5;i++)
						{
							errorLevel="x4.n2 i="+to_string(  i);
							optimumPointsSetMap[i].push_back(generatedPointsLocation[i]);
						}
						errorLevel="x4.n2.1";
						costReciprocalOfOptimum.push_back(1/layoutCost);
						newacceptedCandidate=true;
						errorLevel="x5";
					}else
					{				
						errorLevel="x4.m";
						criterian=cal.newMetrapolisCriterian(iterations,d,costChange,minCostChange,0,T_min,T_max);
						
						if(c==1){
							metroplisCreterian+=to_string(criterian)+"\n";
							iterationsRecord[1]+=to_string(d)+"\n";
							metroplisCreterian2+=to_string(criterian)+"\n";
							if(f==0)
							{
								iterationsRecord[2]+=to_string(d)+"\n";
								metroplisCreterian_300+=to_string(criterian)+"\n";
								metroplisCreterian2_300+=to_string(criterian)+"\n";
							}
						}
						condition=	criterian>cal.randomInRange(0,1); 
						if(condition)
						{
							for(int i=0;i<  5;i++)
							{
								optimumPointsSetMap[i].push_back(generatedPointsLocation[i]);
							}
							costReciprocalOfOptimum.push_back(1/layoutCost);
							newacceptedCandidate=true;
						}
					}

					errorLevel="x4.o";
					if(c==1){
						for (int i = 0; i < 5; i++)
						{
							meansRecord[i]+=to_string(pointMeanPointer[i])+"\n";
							if(f==0)
							{
								meansRecord_300[i]+=to_string( pointMeanPointer[i])+"\n";
							}
						}

					}


					if(newacceptedCandidate)
					{
						if(c==1)
						{
							acceptedCost+=to_string(layoutCost)+"\n";
							iterationsRecord[3]+=to_string(d)+"\n";
						}



						errorLevel="x4.k"+to_string(d);


						//erase elements outside the range of candidate selection pool
						if(costReciprocalOfOptimum.size()>1)
						{
							//messageInfo(cal.matrixToOneLine(costReciprocalOfOptimum));
							if(!isRandomPhase && f%3==0)
							{
								double currentRangeSize=2*(6+k*((((double)domainSize-1)/2)-6));

								for (int ii = 0; ii < costReciprocalOfOptimum.size(); ii++)
								{
									for (int jj = 0; jj < 5; jj++)
									{
										if( (double)optimumPointsSetMap[jj][ii]>= (pointMeanPointer[jj]+(currentRangeSize/2)) || (double)optimumPointsSetMap[jj][ii]<= (pointMeanPointer[jj]-(currentRangeSize/2)))
										{
											/*UpdateResultSummary(to_string(currentRangeSize));
											UpdateResultSummary("d="+to_string(d)+", f="+to_string(f));
											UpdateResultSummary("delete out side range");*/
											for(int gg=0;gg<  5;gg++)
											{
												optimumPointsSetMap[gg].erase(optimumPointsSetMap[gg].begin()+ii);
											}
											costReciprocalOfOptimum.erase(costReciprocalOfOptimum.begin()+ii);
											break;
										}
									}
								}
							}
						}





						//delete the worst element in the group or the oldest
						if(optimumPointsSetMap[0].size()>subGroupSize)
						{	if(f%2==0)
						{
							//get the index of the worst element in the group
							int worstElementIndex = std::min_element(costReciprocalOfOptimum.begin(),costReciprocalOfOptimum.end()) - costReciprocalOfOptimum.begin();
							//UpdateResultSummary("delete worst element");
							for(int i=0;i<  5;i++)
							{
								optimumPointsSetMap[i].erase(optimumPointsSetMap[i].begin()+worstElementIndex);
							}
							costReciprocalOfOptimum.erase(costReciprocalOfOptimum.begin()+worstElementIndex);
						}else
						{
							for(int i=0;i<  5;i++)
							{
								optimumPointsSetMap[i].erase(optimumPointsSetMap[i].begin());
							}
							costReciprocalOfOptimum.erase(costReciprocalOfOptimum.begin());
						}
						}




						sizeOfGuideGroup+=to_string(costReciprocalOfOptimum.size())+"\n";

						//changing the mean
						if(costReciprocalOfOptimum.size()>0)
						{
							//method 1 and 2 (simulatenous calculation)

							vector<double> sum(newPoint.size(),0.0); 
							double sumx=0.0;
							vector<double> pointMeanPointer_method2(5,0);
							vector<double> sum2(newPoint.size(),0.0); 

							for(int i=0;i<optimumPointsSetMap[0].size();i++)
							{
								//Fittness value
								for(int s=0;s<  5;s++)
								{
									sum[s]+=(optimumPointsSetMap[s][i]*costReciprocalOfOptimum[i]);
									sum2[s]+=cal.halfRangeShifting(optimumPointsSetMap[s][i],domainSize-1)*costReciprocalOfOptimum[i];
								}
								sumx+=costReciprocalOfOptimum[i];


							}	
							for(int s=0;s<  5;s++)
							{
								pointMeanPointer[s]=(sum[s]/sumx); 
								pointMeanPointer_method2[s]=cal.halfRangeShifting_back(sum2[s]/sumx,domainSize-1);

							}



							//decide which method to use based on the closest distance
							vector<double> dist1(5,0),dist2(5,0);
							for (int i=0;i<optimumPointsSetMap[0].size();i++)
							{
								for(int s=0;s<  5;s++)
								{
									dist1[s]+=abs(pointMeanPointer[s]-(double)optimumPointsSetMap[s][i])*costReciprocalOfOptimum[i];
									dist2[s]+=abs(cal.halfRangeShifting(pointMeanPointer_method2[s],domainSize-1)-cal.halfRangeShifting(optimumPointsSetMap[s][i],domainSize-1))*costReciprocalOfOptimum[i];
									//make decision when reaching the last point
									if(i==(optimumPointsSetMap[0].size()-1) && dist2[s]<dist1[s])
									{
										pointMeanPointer[s]=pointMeanPointer_method2[s];
									}

								}
							}







						}

					}
				}
				errorLevel="x4.l"+to_string(d);
			}
			if(c==1){
				for (int i = 0; i < 5; i++)
				{
					saveStringToFile(&meansRecord[i],"meansRecord"+to_string(i));
					saveStringToFile(&meansRecord_300[i],"meansRecord_300 "+to_string(i));
					saveStringToFile(&iterationsRecord[i],"iterationsRecord"+to_string(i));
				}
				for (int i = 0; i < 6; i++)
				{
					saveStringToFile(&zscoreRecords[i],"zscoreRecords"+to_string(i));
					saveStringToFile(&zscoreRecords_300[i],"zscoreRecords_300 "+to_string(i));
					saveStringToFile(&subCostsRecords[i],"subCostsRecords"+to_string(i));
				}
				saveStringToFile(&sizeOfGuideGroup,"sizeOfGuideGroup");
				saveStringToFile(&costRecord,"costRecord");
				saveStringToFile(&costRecord_300,"costRecord_300");
				saveStringToFile(&acceptedCost,"acceptedCost");
				saveStringToFile(&metroplisCreterian,"metroplisCreterian");
				saveStringToFile(&metroplisCreterian_300,"metroplisCreterian_300");
				saveStringToFile(&metroplisCreterian2,"metroplisCreterian2");
				saveStringToFile(&metroplisCreterian2_300,"metroplisCreterian2_300");

				UpdateResultSummary("Limda 1 at each load point:");
				UpdateResultSummary(cal.matrixToOneLine(limda1_optimum));
				UpdateResultSummary("Limda 2 at each load point:");
				UpdateResultSummary(cal.matrixToOneLine(limda2_optimum));
				UpdateResultSummary("Limda 3 at each load point:");
				UpdateResultSummary(cal.matrixToOneLine(limda3_optimum));
				UpdateResultSummary("OptimumCosts:");
				UpdateResultSummary(cal.matrixToOneLine(OptimumCosts));
				UpdateResultSummary("OptimumZscore:");
				UpdateResultSummary(cal.matrixToOneLine(OptimumZscore));
				UpdateResultSummary("Indexes of final mean: "+cal.matrixToOneLine(pointMeanPointer));
			}
terminateSA:
			UpdateProgressBar("Optimization finished ");
		}
		return optimumPointsLocation;
	}catch(exception ex)
	{
		result.push_back(errorLevel);
		UpdateResultSummary();
	}
}
void fixturemodule::discretizationUnit(int** &integerOfPointsCloud,vector<	vector<arrayLocation>>&outputPath,vector<arrayLocation> *extremePoints,vector <double>  *weights_x,vector <double> *weights_y)
{
	UpdateProgressBar("Initiate discretization unit");
	faceToPointsCloud(integerOfPointsCloud); 
	int** testCloud;
	copy(testCloud, integerOfPointsCloud);
	UpdateProgressBar("Extract edges");
	isolateBorder(testCloud);
	getBorderPaths(testCloud,outputPath);
	UpdateResultSummary("Number of closed loop paths : "+to_string(outputPath.size()));	
	///////////////////////////////////////////////////////////////////////////////////////////
	UpdateProgressBar("Get reference points");
	getExtremePoints(outputPath[0],*extremePoints,*weights_x,*weights_y);
	numberOfextremePointsOfOuterPath=(*extremePoints).size();
	//maximum distance between two extrem points
	if(!extremPointsAtOuterBoundary)
	{
		(*extremePoints).clear();
		(*weights_x).clear();
		(*weights_y).clear();
		numberOfextremePointsOfOuterPath=0;
	}
	if(outputPath.size()>1)
	{
		for(int i=1;i<outputPath.size();i++)
		{
			vector<arrayLocation> extremePoints_temp;
			vector <double>  weights_x_temp;vector <double> weights_y_temp;
			getExtremePoints(outputPath[i],extremePoints_temp,weights_x_temp,weights_y_temp);
			(*extremePoints).insert((*extremePoints).end(), extremePoints_temp.begin(), extremePoints_temp.end());
			(*weights_x).insert((*weights_x).end(), weights_x_temp.begin(), weights_x_temp.end());
			(*weights_y).insert((*weights_y).end(), weights_y_temp.begin(), weights_y_temp.end());
		}
	}
	//normalize the weights of extreme points
	cal.linearNormilaization(*weights_x);
	cal.linearNormilaization(*weights_y);
}
bool fixturemodule::mathmaticalAttributesUnit(int** &integerOfPointsCloud,vector<	vector<arrayLocation>>&outputPath,vector<arrayLocation> *extremePoints,vector <double>  *weights_x,vector <double> *weights_y,vector <double>&curveture)
{
	string errorLevel;
	try{
		/////////////////////////////////////////////////
		UpdateProgressBar("Compute mathmatical attributes");
		if(k_obj[2]!=0.0)
		{
			distanceScore(integerOfPointsCloud,&outputPath);
			cavitiesScore(*extremePoints, *weights_x,  *weights_y,&outputPath);
		}
		errorLevel="1";
		//////////////////////////////////////////////////
		const	int numberOfPaths =outputPath.size();  
		vector <double>firstDervative;
		vector <double>secondDervative;
		vector <Point2d>backpoint,forwardpoint;
		UpdateProgressBar("Compute mathmatical attributes: Dervatives");
		if(!Dervatives(&outputPath,numberOfPairs,ebsilon,backpoint,forwardpoint,firstDervative,secondDervative,curveture))
		{
			return false;
		}
		if(firstDervative.size()!=outputPath[0].size() ||secondDervative.size()!=outputPath[0].size()||curveture.size()!=outputPath[0].size())
		{
			UpdateResultSummary("Fatal error: attributes are not in the same size of the pont set domain");
			return false;
		}
		UpdateProgressBar("Compute mathmatical attributes: Dervatives");
		errorLevel="2";
		cal.SmoothByBin(curveture,10);
		errorLevel="3";
		//vector <double>firstDervative2_rounded;
		//vector <double>secondDervative2_rounded;
		vector <double>curveture2_rounded;
		cal.SmoothByBin(firstDervative,2);
		cal.SmoothByBin(secondDervative,2);

		errorLevel="4";
		for(int i=0;i<curveture.size();i++)
		{
			//firstDervative2_rounded.push_back( (double)((int)(100*firstDervative[i]))/100);
			//secondDervative2_rounded.push_back( (double)((int)(100*secondDervative[i]))/100);
			curveture2_rounded.push_back( (double)((int)(100*curveture[i]))/100);
		}
		saveStringToFile(&cal.matrixToString( curveture2_rounded),"curveture2_rounded");
		errorLevel="5";
		//firstDervativeEntropy=cal.enyropy(firstDervative2_rounded); messageInfo(firstDervativeEntropy);
		//secondDervativeEntropy=cal.enyropy(secondDervative2_rounded); messageInfo(secondDervativeEntropy);
		curvatureEntropy=cal.enyropy(curveture2_rounded); //messageInfo(curvatureEntropy);
		UpdateResultSummary("curvature entropy = "+to_string(curvatureEntropy));
		iterations=(1+curvatureEntropy*0.4)*iterations;// messageInfo(iterations);
		UpdateResultSummary("Number of iterations: "+to_string( iterations));
		UpdateResultSummary("Diagonal: "+to_string( boundaryDiagonal));
		UpdateProgressBar("Compute mathmatical attributes: Norms");
		errorLevel="6";
		pathPointsUnitNorm_x_all.resize(numberOfPaths);
		pathPointsUnitNorm_y_all.resize(numberOfPaths);
		errorLevel="7";
		for(int i=0;i<outputPath[0].size();i++)
		{
			double xx,yy;
			yy=-(	forwardpoint[i].X-backpoint[i].X);
			xx=(		forwardpoint[i].Y-backpoint[i].Y);
			pathPointsUnitNorm_y_all[0].push_back(yy/sqrt((yy*yy)+(xx*xx)));
			pathPointsUnitNorm_x_all[0].push_back(xx/sqrt((yy*yy)+(xx*xx)));
		}
		errorLevel="8";
		cal.SmoothByBin(pathPointsUnitNorm_y_all[0],2);
		cal.SmoothByBin(pathPointsUnitNorm_x_all[0],2);
		errorLevel="9";
		//Normalize the norms
		for (int i = 0; i < pathPointsUnitNorm_y_all[0].size(); i++)
		{
			double hypotenuse=sqrt(pathPointsUnitNorm_x_all[0][i]*pathPointsUnitNorm_x_all[0][i]+pathPointsUnitNorm_y_all[0][i]*pathPointsUnitNorm_y_all[0][i]);
			pathPointsUnitNorm_x_all[0][i]=pathPointsUnitNorm_x_all[0][i]/hypotenuse;
			pathPointsUnitNorm_y_all[0][i]=pathPointsUnitNorm_y_all[0][i]/hypotenuse;
		}
		UpdateProgressBar("Compute mathmatical attributes: Points intensity");
		if(k_obj[2]!=0.0)
		{
			conicalBeamAndShearLine(integerOfPointsCloud,&outputPath);
		}
		meanDistance=boundaryDiagonal/2;
		return true;
	}catch(exception ex)
	{
		UpdateResultSummary("Error while computing mathmatical attributes");
		UpdateResultSummary("Error level ->"+errorLevel);
		return false;
	}
}
bool fixturemodule::extractSelectionDomain(vector <double>&curveture_all,vector <double>&curveture_subSet,vector<	vector<arrayLocation>>&outputPath)
{
	string errorLevel;
	try{
		//set a subset for the selectionDomain
		UpdateProgressBar("Extract selection domain");
		vector <double>curveture_abs=curveture_all;
		for (int i = 0; i < curveture_abs.size(); i++)
		{
			curveture_abs[i]=abs(curveture_abs[i]);
		}
		std::sort(curveture_abs.rbegin(), curveture_abs.rend());
		errorLevel="1";
		double thresholdCurvatureValue=curveture_abs[(int)(dropPercentage*(double)curveture_abs.size())];
		int dropCounter=0;
		//Check if appropriate beta value
		if((double)beta/(double)outputPath[0].size()>maximumBetaToPointsSize)
		{
			messageInfo("Please specify a smaller beta parameter" );
			return false;
		}
		errorLevel="2";
		for(int i=0;i<outputPath[0].size();i++)
		{
			//protection code
			if(pathPointsUnitNorm_x_all[0][i]!=pathPointsUnitNorm_x_all[0][i])
			{
				UpdateResultSummary("Fatal error: indetermined value detected for path norm_x");
			}
			if(pathPointsUnitNorm_y_all[0][i]!=pathPointsUnitNorm_y_all[0][i])
			{
				UpdateResultSummary("Fatal error: indetermined value detected for path norm_y");
			}
			if(curveture_all[i]!=curveture_all[i])
			{
				UpdateResultSummary("Fatal error: indetermined value detected for path curveture");
			}
			if(outputPath[0].size()!=pathPointsUnitNorm_x_all[0].size() ||
				outputPath[0].size()!=pathPointsUnitNorm_y_all[0].size() ||
				outputPath[0].size()!=curveture_all.size())
			{
				UpdateResultSummary("Fatal error: Mathmatical attributes are not in the same size with point domain");
			}
			errorLevel="3";
			//singularity
			if(abs(curveture_all[i])>thresholdCurvatureValue+ebsilon)
				//if(false)
			{
				errorLevel="4";
				//dropCounter=beta;
				continue;
			}else
			{
				errorLevel="5";
				// spacing
				if(dropCounter==beta)
				{
					dropCounter=0;
				}else
				{
					dropCounter+=1;
					continue;	
				}
				errorLevel="6";
				selectionDomain.push_back(outputPath[0][i]);
				errorLevel="6.1";
				pathPointsUnitNorm_x_subSet.push_back(pathPointsUnitNorm_x_all[0][i]);
				errorLevel="6.2";
				pathPointsUnitNorm_y_subSet.push_back(pathPointsUnitNorm_y_all[0][i]);
				errorLevel="6.3";
				curveture_subSet.push_back(curveture_all[i]);
				errorLevel="6.4";
				if(k_obj[2]!=0.0)
				{
					//These attributes related to the evaluation of the ANN model
					pointsInSmallregion_subSet.push_back(pointsInSmallregion_all[i]);
					errorLevel="6.5";
					pointsInBigRegion_subSet.push_back(pointsInBigRegion_all[i]);
					errorLevel="7";
					minimumDistanceToCavity_subSet.push_back(minimumDistanceToCavity_all[i]);
					cumulativeDistanceToCavity_subSet.push_back(cumulativeDistanceToCavity_all[i]);
					pointsInForawrdConicalBeam_subSet.push_back(pointsInForawrdConicalBeam_all[i]);
					lengthOfShearLine_subSet.push_back(lengthOfShearLine_all[i]);
				}
				errorLevel="8";
			}
		}
		return true;
	}catch(exception ex)
	{
		UpdateResultSummary("Error when extracting the selection domain");
		UpdateResultSummary("Error level ->"+errorLevel);
		return false;
	}
}
