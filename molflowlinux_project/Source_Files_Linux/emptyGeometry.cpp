void LinuxEmptyGeometry() {

	Geometry *geom = worker.GetGeometry();
	mApp->ResetSimulation(false);

	try {
		geom->EmptyGeometry();
		worker.CalcTotalOutgassing();
		//default values
		worker.enableDecay = false;
		worker.halfLife = 1;
		worker.gasMass = 28;
		worker.ResetMoments();
	}
	catch (Error &e) {
		//GLMessageBox::Display((char *)e.GetMsg(), "Error resetting geometry", GLDLG_OK, GLDLG_ICONERROR);
		geom->Clear();
		return;
	}
	//worker.SetCurrentFileName("");
	nbDesStart = 0;
	nbHitStart = 0;

	/*for (int i = 0; i < MAX_VIEWER; i++)
		viewer[i]->SetWorker(&worker);

	//UpdateModelParams();
	startSimu->SetEnabled(true);
	compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(true);
	ClearFacetParams();
	ClearFormulas();
	ClearParameters();
	ClearAllSelections();
	ClearAllViews();

	/*GLParser *f = new GLParser();
	f->SetExpression("A2/SUMDES");
	f->SetName("Trans. Prob.");
	f->Parse();*/
	//AddFormula("Trans.prob.", "A2/SUMDES");

	/*UpdateStructMenu();
	// Send to sub process
	worker.Reload();

	//UpdatePlotters();

	if (timeSettings) timeSettings->RefreshMoments();
	if (momentsEditor) momentsEditor->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (texturePlotter) texturePlotter->Update(0.0, true);
	//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
	if (outgassingMap) outgassingMap->Update(m_fTime, true);
	if (movement) movement->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	if (formulaEditor) formulaEditor->Refresh();
	
	if (textureScaling) textureScaling->Update();
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	
	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
	UpdatePlotters();
	*/
}


/* Original Windows Code
 void MolFlow::EmptyGeometry() {

	Geometry *geom = worker.GetGeometry();
	ResetSimulation(false);

	try {
		geom->EmptyGeometry();
		worker.CalcTotalOutgassing();
		//default values
		worker.enableDecay = false;
		worker.halfLife = 1;
		worker.gasMass = 28;
		worker.ResetMoments();
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error resetting geometry", GLDLG_OK, GLDLG_ICONERROR);
		geom->Clear();
		return;
	}
	worker.SetCurrentFileName("");
	nbDesStart = 0;
	nbHitStart = 0;

	for (int i = 0; i < MAX_VIEWER; i++)
		viewer[i]->SetWorker(&worker);

	//UpdateModelParams();
	startSimu->SetEnabled(true);
	compACBtn->SetEnabled(modeCombo->GetSelectedIndex() == 1);
	//resetSimu->SetEnabled(true);
	ClearFacetParams();
	ClearFormulas();
	ClearParameters();
	ClearAllSelections();
	ClearAllViews();

	/*GLParser *f = new GLParser();
	f->SetExpression("A2/SUMDES");
	f->SetName("Trans. Prob.");
	f->Parse();*/
	//AddFormula("Trans.prob.", "A2/SUMDES");

	/* UpdateStructMenu();
	// Send to sub process
	worker.Reload();

	//UpdatePlotters();

	if (timeSettings) timeSettings->RefreshMoments();
	if (momentsEditor) momentsEditor->Refresh();
	if (pressureEvolution) pressureEvolution->Refresh();
	if (timewisePlotter) timewisePlotter->Refresh();
	if (profilePlotter) profilePlotter->Refresh();
	if (texturePlotter) texturePlotter->Update(0.0, true);
	//if (parameterEditor) parameterEditor->UpdateCombo(); //Done by ClearParameters()
	if (outgassingMap) outgassingMap->Update(m_fTime, true);
	if (movement) movement->Update();
	if (globalSettings && globalSettings->IsVisible()) globalSettings->Update();
	if (formulaEditor) formulaEditor->Refresh();
	
	if (textureScaling) textureScaling->Update();
	if (facetDetails) facetDetails->Update();
	if (facetCoordinates) facetCoordinates->UpdateFromSelection();
	if (vertexCoordinates) vertexCoordinates->Update();
	
	UpdateTitle();
	changedSinceSave = false;
	ResetAutoSaveTimer();
	UpdatePlotters();
}
 */
