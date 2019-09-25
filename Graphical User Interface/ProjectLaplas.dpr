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
  Unitamg1r5Parameters in 'Unitamg1r5Parameters.pas' {Formamg1r5Parameters};

{$R *.res}

begin
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
  Application.Run;
end.
