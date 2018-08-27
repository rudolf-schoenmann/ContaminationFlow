void LinuxLoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("molflow.cfg");
		MolflowGeometry *geom = worker.GetMolflowGeometry();

		/*
        f->ReadKeyword("showRules"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("showTimeOverlay"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTime = f->ReadInt();
		f->ReadKeyword("texColormap"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			//viewer[i]->showColormap = 
			f->ReadInt();
		f->ReadKeyword("translation"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumHits = f->ReadLLong();
		f->ReadKeyword("dispNumLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumLeaks = f->ReadLLong();
		f->ReadKeyword("dirShow"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();*/
		f->ReadKeyword("dirNorme"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)f->ReadDouble());
		f->ReadKeyword("dirAutoNormalize"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("dirCenter"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		/*f->ReadKeyword("angle"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->angleStep = f->ReadDouble();*/
		f->ReadKeyword("autoScale"); f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("autoScale_include_constant_flow"); f->ReadKeyword(":");
		geom->texAutoScaleIncludeConstantFlow = f->ReadInt();

		f->ReadKeyword("textures_min_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("processNum"); f->ReadKeyword(":");
		nbProc = f->ReadLLong();
#ifdef _DEBUG
		nbProc = 1;
#endif
		if (nbProc <= 0) nbProc = 1;
		f->ReadKeyword("recents"); f->ReadKeyword(":"); f->ReadKeyword("{");
		w = f->ReadString();
		while (strcmp(w, "}") != 0 && nbRecent < MAX_RECENT) {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}

		f->ReadKeyword("cdir"); f->ReadKeyword(":");
		strcpy(currentDir, f->ReadString());
		f->ReadKeyword("cseldir"); f->ReadKeyword(":");
		strcpy(currentSelDir, f->ReadString());
		f->ReadKeyword("autonorme"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("autoSaveFrequency"); f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly"); f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates"); f->ReadKeyword(":");
		/*checkForUpdates =*/ f->ReadInt(); //Old checkforupdates
		f->ReadKeyword("autoUpdateFormulas"); f->ReadKeyword(":");
		autoUpdateFormulas = f->ReadInt();
		f->ReadKeyword("compressSavedFiles"); f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();
		f->ReadKeyword("gasMass"); f->ReadKeyword(":");
		worker.gasMass = f->ReadDouble();
		f->ReadKeyword("expandShortcutPanel"); f->ReadKeyword(":");
		bool isOpen = f->ReadInt();
		if (isOpen) shortcutPanel->Open();
		else shortcutPanel->Close();
		f->ReadKeyword("hideLot"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->hideLot = f->ReadInt();
		f->ReadKeyword("lowFluxMode"); f->ReadKeyword(":");
		worker.ontheflyParams.lowFluxMode = f->ReadInt();
		f->ReadKeyword("lowFluxCutoff"); f->ReadKeyword(":");
		worker.ontheflyParams.lowFluxCutoff = f->ReadDouble();
		f->ReadKeyword("leftHandedView"); f->ReadKeyword(":");
		leftHandedView = f->ReadInt();
	}
	catch (...) {
		/*std::ostringstream tmp;
		tmp << err.GetMsg() << "\n\nThis is normal on the first launch and if you upgrade from an earlier version\n";
		tmp << "MolFlow will use default program settings.\nWhen you quit, a correct config file will be written\n";
		GLMessageBox::Display(tmp.str().c_str(), "Error loading config file", GLDLG_OK, GLDLG_ICONINFO);*/
		SAFE_DELETE(f);
		return;
	}
	SAFE_DELETE(f);
}


/* Old Windows Code to recycle
void MolFlow::LoadConfig() {

	FileReader *f = NULL;
	char *w;
	nbRecent = 0;

	try {

		f = new FileReader("molflow.cfg");
		MolflowGeometry *geom = worker.GetMolflowGeometry();

		f->ReadKeyword("showRules"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showRule = f->ReadInt();
		f->ReadKeyword("showNormals"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showNormal = f->ReadInt();
		f->ReadKeyword("showUV"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showUV = f->ReadInt();
		f->ReadKeyword("showLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLine = f->ReadInt();
		f->ReadKeyword("showLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showLeak = f->ReadInt();
		f->ReadKeyword("showHits"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHit = f->ReadInt();
		f->ReadKeyword("showVolume"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVolume = f->ReadInt();
		f->ReadKeyword("showTexture"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTexture = f->ReadInt();
		f->ReadKeyword("showFilter"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showFilter = f->ReadInt();
		f->ReadKeyword("showIndices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showIndex = f->ReadInt();
		f->ReadKeyword("showVertices"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showVertex = f->ReadInt();
		f->ReadKeyword("showMode"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showBack = f->ReadInt();
		f->ReadKeyword("showMesh"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showMesh = f->ReadInt();
		f->ReadKeyword("showHidden"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHidden = f->ReadInt();
		f->ReadKeyword("showHiddenVertex"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showHiddenVertex = f->ReadInt();
		f->ReadKeyword("showTimeOverlay"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showTime = f->ReadInt();
		f->ReadKeyword("texColormap"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			//viewer[i]->showColormap = 
			f->ReadInt();
		f->ReadKeyword("translation"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->transStep = f->ReadDouble();
		f->ReadKeyword("dispNumLines"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumHits = f->ReadLLong();
		f->ReadKeyword("dispNumLeaks"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->dispNumLeaks = f->ReadLLong();
		f->ReadKeyword("dirShow"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->showDir = f->ReadInt();
		f->ReadKeyword("dirNorme"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)f->ReadDouble());
		f->ReadKeyword("dirAutoNormalize"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("dirCenter"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("angle"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->angleStep = f->ReadDouble();
		f->ReadKeyword("autoScale"); f->ReadKeyword(":");
		geom->texAutoScale = f->ReadInt();
		f->ReadKeyword("autoScale_include_constant_flow"); f->ReadKeyword(":");
		geom->texAutoScaleIncludeConstantFlow = f->ReadInt();

		f->ReadKeyword("textures_min_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_all"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_pressure_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[0].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_all"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_impingement_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[1].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("textures_min_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.all = f->ReadDouble();
		f->ReadKeyword("textures_min_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.min.moments_only = f->ReadDouble();
		f->ReadKeyword("textures_max_density_all"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.all = f->ReadDouble();
		f->ReadKeyword("textures_max_density_moments_only"); f->ReadKeyword(":");
		geom->texture_limits[2].autoscale.max.moments_only = f->ReadDouble();

		f->ReadKeyword("processNum"); f->ReadKeyword(":");
		nbProc = f->ReadLLong();
#ifdef _DEBUG
		nbProc = 1;
#endif
		if (nbProc <= 0) nbProc = 1;
		f->ReadKeyword("recents"); f->ReadKeyword(":"); f->ReadKeyword("{");
		w = f->ReadString();
		while (strcmp(w, "}") != 0 && nbRecent < MAX_RECENT) {
			recents[nbRecent] = _strdup(w);
			nbRecent++;
			w = f->ReadString();
		}

		f->ReadKeyword("cdir"); f->ReadKeyword(":");
		strcpy(currentDir, f->ReadString());
		f->ReadKeyword("cseldir"); f->ReadKeyword(":");
		strcpy(currentSelDir, f->ReadString());
		f->ReadKeyword("autonorme"); f->ReadKeyword(":");
		geom->SetAutoNorme(f->ReadInt());
		f->ReadKeyword("centernorme"); f->ReadKeyword(":");
		geom->SetCenterNorme(f->ReadInt());
		f->ReadKeyword("normeratio"); f->ReadKeyword(":");
		geom->SetNormeRatio((float)(f->ReadDouble()));
		f->ReadKeyword("autoSaveFrequency"); f->ReadKeyword(":");
		autoSaveFrequency = f->ReadDouble();
		f->ReadKeyword("autoSaveSimuOnly"); f->ReadKeyword(":");
		autoSaveSimuOnly = f->ReadInt();
		f->ReadKeyword("checkForUpdates"); f->ReadKeyword(":");
		/*checkForUpdates =*/ /*f->ReadInt(); //Old checkforupdates
		f->ReadKeyword("autoUpdateFormulas"); f->ReadKeyword(":");
		autoUpdateFormulas = f->ReadInt();
		f->ReadKeyword("compressSavedFiles"); f->ReadKeyword(":");
		compressSavedFiles = f->ReadInt();
		f->ReadKeyword("gasMass"); f->ReadKeyword(":");
		worker.gasMass = f->ReadDouble();
		f->ReadKeyword("expandShortcutPanel"); f->ReadKeyword(":");
		bool isOpen = f->ReadInt();
		if (isOpen) shortcutPanel->Open();
		else shortcutPanel->Close();
		f->ReadKeyword("hideLot"); f->ReadKeyword(":");
		for (int i = 0; i < MAX_VIEWER; i++)
			viewer[i]->hideLot = f->ReadInt();
		f->ReadKeyword("lowFluxMode"); f->ReadKeyword(":");
		worker.ontheflyParams.lowFluxMode = f->ReadInt();
		f->ReadKeyword("lowFluxCutoff"); f->ReadKeyword(":");
		worker.ontheflyParams.lowFluxCutoff = f->ReadDouble();
		f->ReadKeyword("leftHandedView"); f->ReadKeyword(":");
		leftHandedView = f->ReadInt();
	}
	catch (...) {
		/*std::ostringstream tmp;
		tmp << err.GetMsg() << "\n\nThis is normal on the first launch and if you upgrade from an earlier version\n";
		tmp << "MolFlow will use default program settings.\nWhen you quit, a correct config file will be written\n";
		GLMessageBox::Display(tmp.str().c_str(), "Error loading config file", GLDLG_OK, GLDLG_ICONINFO);*/
		SAFE_DELETE(f);
		return;
	}
	SAFE_DELETE(f);
}*/
