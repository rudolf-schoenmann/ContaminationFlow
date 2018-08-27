void LinuxLoadFile(char *fName) {	
	/*
    char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open file", fileLFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}

	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);*/
    
    char fullName[512];
    strcpy(fullName, fName);
    
    
	try {
		/*
        ClearFormulas();
		ClearParameters();
		ClearAllSelections();
		ClearAllViews();
		ResetSimulation(false);
		*/

		worker.LoadGeometry(fullName);

		Geometry *geom = worker.GetGeometry();

		
		// Default initialisation
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->SetWorker(&worker);

		//UpdateModelParams();
		startSimu->SetEnabled(true);
		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		ClearFacetParams();
		nbDesStart = worker.nbDesorption;
		nbHitStart = worker.nbMCHit;
		AddRecent(fullName);
		geom->viewStruct = -1;

		UpdateStructMenu();
		UpdateCurrentDir(fullName);

		// Check non simple polygon
		progressDlg2->SetMessage("Checking for non simple polygons...");

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();
		// Set up view
		// Default
		/*
        viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);
		

		ResetAutoSaveTimer(); 
		//UpdatePlotters();

		if (timeSettings) timeSettings->RefreshMoments();
		if (momentsEditor) momentsEditor->Refresh();
		if (pressureEvolution) pressureEvolution->Refresh();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (profilePlotter) profilePlotter->Refresh();
		if (texturePlotter) texturePlotter->Update(0.0,true);
		//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
		if (textureScaling) textureScaling->Update();
		if (outgassingMap) outgassingMap->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (movement) movement->Update();
		if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
		if (formulaEditor) formulaEditor->Refresh();
		*/
	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		/*GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);*/
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(false); //?
	SAFE_DELETE(progressDlg2); //?
	changedSinceSave = false; //?
}




/* Old Windows Code to recycle
 * 
void MolFlow::LoadFile(char *fName) {

	char fullName[512];
	char shortName[512];
	strcpy(fullName, "");

	if (fName == NULL) {
		FILENAME *fn = GLFileBox::OpenFile(currentDir, NULL, "Open file", fileLFilters, 0);
		if (fn)
			strcpy(fullName, fn->fullName);
	}
	else {
		strcpy(fullName, fName);
	}

	GLProgress *progressDlg2 = new GLProgress("Preparing to load file...", "Please wait");
	progressDlg2->SetVisible(true);
	progressDlg2->SetProgress(0.0);
	//GLWindowManager::Repaint();

	if (strlen(fullName) == 0) {
		progressDlg2->SetVisible(false);
		SAFE_DELETE(progressDlg2);
		return;
	}

	char *lPart = strrchr(fullName, '\\');
	if (lPart) strcpy(shortName, lPart + 1);
	else strcpy(shortName, fullName);

	try {
		ClearFormulas();
		ClearParameters();
		ClearAllSelections();
		ClearAllViews();
		ResetSimulation(false);

		worker.LoadGeometry(fullName);

		Geometry *geom = worker.GetGeometry();

		
		// Default initialisation
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->SetWorker(&worker);

		//UpdateModelParams();
		startSimu->SetEnabled(true);
		compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		singleACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
		//resetSimu->SetEnabled(true);
		ClearFacetParams();
		nbDesStart = worker.nbDesorption;
		nbHitStart = worker.nbMCHit;
		AddRecent(fullName);
		geom->viewStruct = -1;

		UpdateStructMenu();
		UpdateCurrentDir(fullName);

		// Check non simple polygon
		progressDlg2->SetMessage("Checking for non simple polygons...");

		geom->CheckCollinear();
		geom->CheckNonSimple();
		geom->CheckIsolatedVertex();
		// Set up view
		// Default
		viewer[0]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[0]->ToFrontView();
		viewer[1]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[1]->ToTopView();
		viewer[2]->SetProjection(ORTHOGRAPHIC_PROJ);
		viewer[2]->ToSideView();
		viewer[3]->SetProjection(PERSPECTIVE_PROJ);
		viewer[3]->ToFrontView();
		SelectViewer(0);

		ResetAutoSaveTimer();
		//UpdatePlotters();

		if (timeSettings) timeSettings->RefreshMoments();
		if (momentsEditor) momentsEditor->Refresh();
		if (pressureEvolution) pressureEvolution->Refresh();
		if (timewisePlotter) timewisePlotter->Refresh();
		if (profilePlotter) profilePlotter->Refresh();
		if (texturePlotter) texturePlotter->Update(0.0,true);
		//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
		if (textureScaling) textureScaling->Update();
		if (outgassingMap) outgassingMap->Update(m_fTime, true);
		if (facetDetails) facetDetails->Update();
		if (facetCoordinates) facetCoordinates->UpdateFromSelection();
		if (vertexCoordinates) vertexCoordinates->Update();
		if (movement) movement->Update();
		if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
		if (formulaEditor) formulaEditor->Refresh();
	}
	catch (Error &e) {

		char errMsg[512];
		sprintf(errMsg, "%s\nFile:%s", e.GetMsg(), shortName);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		RemoveRecent(fName);

	}
	progressDlg2->SetVisible(false);
	SAFE_DELETE(progressDlg2);
	changedSinceSave = false;
} */


