program ProjectLaplas;

uses
  Forms,
  VisualUnit in 'VisualUnit.pas' {Laplas},
  cabinetUnit in 'cabinetUnit.pas' {CabinetForm},
  addBlockUnit in 'addBlockUnit.pas' {AddBlockForm},
  AddSourceUnit in 'AddSourceUnit.pas' {AddSourceForm},
  AddWallUnit in 'AddWallUnit.pas' {AddWallForm},
  UnitCopy in 'UnitCopy.pas' {FormCopyObject},
  UnitFormUnion in 'UnitFormUnion.pas' {FormUnion},
  MeshUnit in 'MeshUnit.pas' {MeshForm},
  UnitGravity in 'UnitGravity.pas' {FormGravity},
  UnitSolidMatLib in 'UnitSolidMatLib.pas' {FormSolidLibMat},
  UnitUserDefinedSolidMaterial in 'UnitUserDefinedSolidMaterial.pas' {FormUserDefinedSolidMat},
  UnitFlyuidMatLib in 'UnitFlyuidMatLib.pas' {FormFluidLibMat},
  UnitUserDefinedFluidMaterial in 'UnitUserDefinedFluidMaterial.pas' {FormUserDefinedFluidMaterial},
  UnitEQGD in 'UnitEQGD.pas' {EGDForm},
  UnitSelProjMat in 'UnitSelProjMat.pas' {FormSelProjMat},
  UnitnonNewtonianFluid in 'UnitnonNewtonianFluid.pas' {FormnonNewtonFluid},
  UnitVariables in 'UnitVariables.pas' {FormVariables},
  UnitSmagorinsky in 'UnitSmagorinsky.pas' {FormSmagorinsky},
  UnitPowerList in 'UnitPowerList.pas' {FormPowerList},
  Unitaddinunion in 'Unitaddinunion.pas' {Form_formirateunion},
  AddVariableUnit in 'AddVariableUnit.pas' {AddVariableForm},
  UnitExportSYMMIC in 'UnitExportSYMMIC.pas' {Formexport2SYMMIC},
  Unitrectangularplot in 'Unitrectangularplot.pas' {frmRectangularPlot},
  Unitresidualplot in 'Unitresidualplot.pas' {Formresidual},
  UnitParallelSetting in 'UnitParallelSetting.pas' {FormSetting},
  UnitTransientMenu in 'UnitTransientMenu.pas' {FormUnsteady},
  UnitSquareWave in 'UnitSquareWave.pas' {FormSquareWave},
  UnitSquareWaveAPPARAT in 'UnitSquareWaveAPPARAT.pas' {FormAPPARAT_Square_Wave},
  UnitTimedependpowerLaw in 'UnitTimedependpowerLaw.pas' {FormTransientPowerSetting},
  Unitresidual2 in 'Unitresidual2.pas' {Formresidual2},
  Unitwallinitposition in 'Unitwallinitposition.pas' {Formwallgeometryposition_init},
  Unitdefineindentcabinet in 'Unitdefineindentcabinet.pas' {Formcabinetindent},
  UnitBlockViewFactors in 'UnitBlockViewFactors.pas' {FormRadiation},
  Unitamgmanager in 'Unitamgmanager.pas' {Form_amg_manager},
  Unitusertempdepend in 'Unitusertempdepend.pas' {Formusertempdepend},
  UnitInitialization in 'UnitInitialization.pas' {FormSpeedInitialization},
  Unit_hot_cold in 'Unit_hot_cold.pas' {Form_hot_cold},
  UnitXYPlot in 'UnitXYPlot.pas' {FormXYPlot},
  UnitRenameVariable in 'UnitRenameVariable.pas' {FormRenameVar},
  UnitScale in 'UnitScale.pas' {FormScale},
  UnitAMGCLManager in 'UnitAMGCLManager.pas' {FormAMGCLParameters},
  Unitamg1r5Parameters in 'Unitamg1r5Parameters.pas' {Formamg1r5Parameters},
  UnitresidualPlotSpallartAllmares in 'UnitresidualPlotSpallartAllmares.pas' {FormResidualSpallart_Allmares},
  UnitResidualSATemp2 in 'UnitResidualSATemp2.pas' {FormResidualSATemp},
  UnitResidualMenterSST in 'UnitResidualMenterSST.pas' {FormResidualSST},
  UnitResidualSSTTemperature in 'UnitResidualSSTTemperature.pas' {FormResidualSSTTemp},
  UnitResidualStandartKEpsilon in 'UnitResidualStandartKEpsilon.pas' {FormResidualStandartKEpsilon},
  UnitResidualStandartK_Epsilon_TEMP in 'UnitResidualStandartK_Epsilon_TEMP.pas' {FormResidualStandart_k_epsilon_Temp},
  Unitpiecewiseconst in 'Unitpiecewiseconst.pas' {Formpiecewiseconstant},
  UnitPatternDelete in 'UnitPatternDelete.pas' {FormPattern},
  UnitTextNameSourcePattern in 'UnitTextNameSourcePattern.pas' {FormTextNameSourcePattern},
  UnitSplash in 'UnitSplash.pas' {SplashForm},
  UnitViewFactors in 'UnitViewFactors.pas' {FormViewFactors},
  Unit_debug in 'Unit_debug.pas' {Form_debug_panel},
  UnitPlaneSelect in 'UnitPlaneSelect.pas' {FormSelectPlaneRotation},
  UnitResidual_Langtry_Menter in 'UnitResidual_Langtry_Menter.pas' {FormResidual_Langtry_Menter},
  UnitResidual_Langtry_Menter_Temp in 'UnitResidual_Langtry_Menter_Temp.pas' {FormResidual_Lagtry_Menter_Temp},
  UnitOptimetric in 'UnitOptimetric.pas' {FormOptimetric};

{$R *.res}

begin
   SplashForm:=TSplashForm.Create(nil);
   SplashForm.Show;
   SplashForm.Repaint;


  Application.Initialize;
  // главная форма визуализатор
  Application.CreateForm(TLaplas, Laplas);
  Application.CreateForm(TCabinetForm, CabinetForm);
  Application.CreateForm(TAddBlockForm, AddBlockForm);
  Application.CreateForm(TAddSourceForm, AddSourceForm);
  Application.CreateForm(TAddWallForm, AddWallForm);
  Application.CreateForm(TFormCopyObject, FormCopyObject);
  Application.CreateForm(TFormUnion, FormUnion);
  Application.CreateForm(TMeshForm, MeshForm);
  Application.CreateForm(TFormGravity, FormGravity);
  Application.CreateForm(TFormSolidLibMat, FormSolidLibMat);
  Application.CreateForm(TFormUserDefinedSolidMat, FormUserDefinedSolidMat);
  Application.CreateForm(TFormFluidLibMat, FormFluidLibMat);
  Application.CreateForm(TFormUserDefinedFluidMaterial, FormUserDefinedFluidMaterial);
  Application.CreateForm(TEGDForm, EGDForm);
  Application.CreateForm(TFormSelProjMat, FormSelProjMat);
  Application.CreateForm(TFormnonNewtonFluid, FormnonNewtonFluid);
  Application.CreateForm(TFormVariables, FormVariables);
  Application.CreateForm(TFormSmagorinsky, FormSmagorinsky);
  Application.CreateForm(TFormPowerList, FormPowerList);
  Application.CreateForm(TForm_formirateunion, Form_formirateunion);
  Application.CreateForm(TAddVariableForm, AddVariableForm);
  Application.CreateForm(TFormexport2SYMMIC, Formexport2SYMMIC);
  Application.CreateForm(TfrmRectangularPlot, frmRectangularPlot);
  Application.CreateForm(TFormresidual, Formresidual);
  Application.CreateForm(TFormSetting, FormSetting);
  Application.CreateForm(TFormUnsteady, FormUnsteady);
  Application.CreateForm(TFormSquareWave, FormSquareWave);
  Application.CreateForm(TFormAPPARAT_Square_Wave, FormAPPARAT_Square_Wave);
  Application.CreateForm(TFormTransientPowerSetting, FormTransientPowerSetting);
  Application.CreateForm(TFormresidual2, Formresidual2);
  Application.CreateForm(TFormwallgeometryposition_init, Formwallgeometryposition_init);
  Application.CreateForm(TFormcabinetindent, Formcabinetindent);
  Application.CreateForm(TFormRadiation, FormRadiation);
  Application.CreateForm(TForm_amg_manager, Form_amg_manager);
  Application.CreateForm(TFormusertempdepend, Formusertempdepend);
  Application.CreateForm(TFormSpeedInitialization, FormSpeedInitialization);
  Application.CreateForm(TForm_hot_cold, Form_hot_cold);
  Application.CreateForm(TFormXYPlot, FormXYPlot);
  Application.CreateForm(TFormRenameVar, FormRenameVar);
  Application.CreateForm(TFormScale, FormScale);
  Application.CreateForm(TFormAMGCLParameters, FormAMGCLParameters);
  Application.CreateForm(TFormamg1r5Parameters, Formamg1r5Parameters);
  Application.CreateForm(TFormResidualSpallart_Allmares, FormResidualSpallart_Allmares);
  Application.CreateForm(TFormResidualSATemp, FormResidualSATemp);
  Application.CreateForm(TFormResidualSST, FormResidualSST);
  Application.CreateForm(TFormResidualSSTTemp, FormResidualSSTTemp);
  Application.CreateForm(TFormResidualStandartKEpsilon, FormResidualStandartKEpsilon);
  Application.CreateForm(TFormResidualStandart_k_epsilon_Temp, FormResidualStandart_k_epsilon_Temp);
  Application.CreateForm(TFormpiecewiseconstant, Formpiecewiseconstant);
  Application.CreateForm(TFormPattern, FormPattern);
  Application.CreateForm(TFormTextNameSourcePattern, FormTextNameSourcePattern);
  Application.CreateForm(TFormViewFactors, FormViewFactors);
  Application.CreateForm(TForm_debug_panel, Form_debug_panel);
  Application.CreateForm(TFormSelectPlaneRotation, FormSelectPlaneRotation);
  Application.CreateForm(TFormResidual_Langtry_Menter, FormResidual_Langtry_Menter);
  Application.CreateForm(TFormResidual_Lagtry_Menter_Temp, FormResidual_Lagtry_Menter_Temp);
  Application.CreateForm(TFormOptimetric, FormOptimetric);
  //Application.CreateForm(TSplashForm, SplashForm);

  SplashForm.Hide;
  SplashForm.Free;

  Application.Run;
end.
