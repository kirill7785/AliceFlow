unit UnitParallelSetting;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls, Vcl.AppEvnts;

type
  TFormSetting = class(TForm)
    rgParallel: TRadioGroup;
    PanelSolverSetting: TPanel;
    Labelalgocoupling: TLabel;
    ComboBoxPressureVelocityCoupling: TComboBox;
    Label1: TLabel;
    ComboBoxFlowScheme: TComboBox;
    PanelSolverSetting2: TPanel;
    LabelSolverselect: TLabel;
    ComboBoxSolverSetting: TComboBox;
    Button_amg_manager: TButton;
    ApplicationEvents1: TApplicationEvents;
    Label2: TLabel;
    ComboBoxSchemeTemperature: TComboBox;
    GroupBox_lfil: TGroupBox;
    ComboBox_lfil: TComboBox;
    Label3: TLabel;
    GroupBox1: TGroupBox;
    Label4: TLabel;
    ComboBox_m_restart_for_gmres: TComboBox;
    GroupBox2: TGroupBox;
    ComboBoxStaticStructuralSolverSetting: TComboBox;
    GroupBoxNumberProcessors: TGroupBox;
    ComboBoxNumberProcessors: TComboBox;
    ComboBoxFlowSchemePrefix: TComboBox;
    ComboBoxGPU_id: TComboBox;
    Label5: TLabel;
    GroupBoxChebyshevPolynomdegree: TGroupBox;
    ComboBoxChebyshevdegree: TComboBox;
    LabelChebyshevdegree: TLabel;
    Label6: TLabel;
    ComboBoxSchemeTurbulent: TComboBox;
    // Галочка NonLinearmultipass увеличивает время счёта в 14раз, а
    // точность на модели tgf01 как была так и осталась не меняется с точностью
    // до долей градуса. 11.02.2023 Рекомендую снять галочку на моделях tgf
    // даже с нелинейными свойствами материалов при нестационарном расчёте.
    CheckBoxNonLinearmultipass: TCheckBox;
    procedure Button_amg_managerClick(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure ComboBoxSolverSettingChange(Sender: TObject);
    procedure ComboBoxStaticStructuralSolverSettingChange(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormSetting: TFormSetting;

implementation

{$R *.dfm}

uses Unitamgmanager, MeshUnit, VisualUnit, UnitAMGCLManager,
  Unitamg1r5Parameters, UnitEQGD;


// Запрет форме сворачиваться.
procedure TFormSetting.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
   if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

// Вызывает окно настроек параметров.
procedure TFormSetting.Button_amg_managerClick(Sender: TObject);
begin
    // Метод контрольного объёма
    if ((Laplas.egddata.itemper=1)
    or (Laplas.egddata.myflmod[0].iflow=1)) then
    begin
       if (ComboBoxSolverSetting.ItemIndex=5) then
       begin
          // GPU AMGCL Denis Demidov BiCGStab + samg
          FormAMGCLParameters.ShowModal;
       end
       else if ((ComboBoxSolverSetting.ItemIndex=4)or
                (ComboBoxSolverSetting.ItemIndex=12)or
                (ComboBoxSolverSetting.ItemIndex=13)or
                (ComboBoxSolverSetting.ItemIndex=16)or
                (ComboBoxSolverSetting.ItemIndex=17)) then
       begin
          // 12, 13, 16 and 17 Chebyshev + AMGCL.

          // CPU AMGCL Denis Demidov BiCGStab + samg
          FormAMGCLParameters.ShowModal;
       end
       else  if ((ComboBoxSolverSetting.ItemIndex=1)or
                 (ComboBoxSolverSetting.ItemIndex=11)) then
       begin
          // 11 Chebyshev + amg1r5.

          // amg1r5 Ruge and Stueben.
          Formamg1r5Parameters.ShowModal;
       end
        else  if ((ComboBoxSolverSetting.ItemIndex=14)or
        (ComboBoxSolverSetting.ItemIndex=15)or
        (ComboBoxSolverSetting.ItemIndex=3)) then
       begin
          // Вызов настроек РУМБА v.0.14 для поправки давления.
          // Гибридный алгоритм. Механика неактивна, для скорости BiCGStab+ILU2 + для поправки давления РУМБА v.0.14

          if (Laplas.egddata.itemper=0) then
          begin
             // Уравнение теплопередачи не решается.
             Form_amg_manager.PanelTemperature1.Visible:=false;
             Form_amg_manager.PanelTemperature2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogTemperature.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Temperature.Visible:=false;
             Form_amg_manager.CheckBoxTemperatureMatrixPortrait.Visible:=false;
          end
          else
          begin
             // Уравнение теплопередачи решается.
             Form_amg_manager.PanelTemperature1.Visible:=true;
             Form_amg_manager.PanelTemperature2.Visible:=true;
             Form_amg_manager.CheckBoxprintlogTemperature.Visible:=true;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Temperature.Visible:=true;
             Form_amg_manager.CheckBoxTemperatureMatrixPortrait.Visible:=true;
          end;
          if (Laplas.egddata.myflmod[0].iflow=1) then
          begin
             // cfd is active
             Form_amg_manager.PanelSpeed1.Visible:=false;
             Form_amg_manager.PanelPressure1.Visible:=true;
             Form_amg_manager.PanelSpeed2.Visible:=false;
             Form_amg_manager.PanelPressure2.Visible:=true;
             Form_amg_manager.CheckBoxprintlogSpeed.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=false;
             Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=false;
             Form_amg_manager.CheckBoxprintlogPressure.Visible:=true;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=true;
             Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=true;
          end
            else
          begin
             // cfd not active
             // Поправка давления активна а скорость неактивна.
             Form_amg_manager.PanelSpeed1.Visible:=false;
             Form_amg_manager.PanelPressure1.Visible:=false;
             Form_amg_manager.PanelSpeed2.Visible:=false;
             Form_amg_manager.PanelPressure2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogSpeed.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=false;
             Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=false;
             Form_amg_manager.CheckBoxprintlogPressure.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=false;
             Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=false;
          end;


           // Mechanical not active.
           Form_amg_manager.PanelStress1.Visible:=false;
           Form_amg_manager.PanelStress2.Visible:=false;
           Form_amg_manager.CheckBoxprintlogStress.Visible:=false;
           Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Stress.Visible:=false;
           Form_amg_manager.CheckBoxStressMatrixPortrait.Visible:=false;

           Form_amg_manager.ShowModal;

       end
       else if (ComboBoxSolverSetting.ItemIndex=7) then
       begin
          // Вызов настроек РУМБА 0.14

          if (Laplas.egddata.itemper=0) then
          begin
             // Уравнение теплопередачи не решается.
             Form_amg_manager.PanelTemperature1.Visible:=false;
             Form_amg_manager.PanelTemperature2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogTemperature.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Temperature.Visible:=false;
             Form_amg_manager.CheckBoxTemperatureMatrixPortrait.Visible:=false;
          end
          else
          begin
             // Уравнение теплопередачи решается.
             Form_amg_manager.PanelTemperature1.Visible:=true;
             Form_amg_manager.PanelTemperature2.Visible:=true;
             Form_amg_manager.CheckBoxprintlogTemperature.Visible:=true;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Temperature.Visible:=true;
             Form_amg_manager.CheckBoxTemperatureMatrixPortrait.Visible:=true;
          end;

          if (Laplas.egddata.myflmod[0].iflow=1) then
          begin
             // cfd is active
             Form_amg_manager.PanelSpeed1.Visible:=true;
             Form_amg_manager.PanelPressure1.Visible:=true;
             Form_amg_manager.PanelSpeed2.Visible:=true;
             Form_amg_manager.PanelPressure2.Visible:=true;
             Form_amg_manager.CheckBoxprintlogSpeed.Visible:=true;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=true;
             Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=true;
             Form_amg_manager.CheckBoxprintlogPressure.Visible:=true;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=true;
             Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=true;
          end
            else
          begin
             // cfd not active
             Form_amg_manager.PanelSpeed1.Visible:=false;
             Form_amg_manager.PanelPressure1.Visible:=false;
             Form_amg_manager.PanelSpeed2.Visible:=false;
             Form_amg_manager.PanelPressure2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogSpeed.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=false;
             Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=false;
             Form_amg_manager.CheckBoxprintlogPressure.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=false;
             Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=false;
          end;
          if (Laplas.egddata.iStaticStructural=0) then
          begin
             // Mechanical not active.
             Form_amg_manager.PanelStress1.Visible:=false;
             Form_amg_manager.PanelStress2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogStress.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Stress.Visible:=false;
             Form_amg_manager.CheckBoxStressMatrixPortrait.Visible:=false;

          end;
          Form_amg_manager.ShowModal;
       end;
    end;
    // Метод конечных элементов.
    if (Laplas.egddata.itemper=2) then
    begin
        if (ComboBoxStaticStructuralSolverSetting.ItemIndex=2) then
        begin
           // Вызов настроек РУМБА 0.14
           // cfd not active
           Form_amg_manager.PanelSpeed1.Visible:=false;
           Form_amg_manager.PanelPressure1.Visible:=false;
           Form_amg_manager.PanelSpeed2.Visible:=false;
           Form_amg_manager.PanelPressure2.Visible:=false;
           Form_amg_manager.CheckBoxprintlogSpeed.Visible:=false;
           Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=false;
           Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=false;
           Form_amg_manager.CheckBoxprintlogPressure.Visible:=false;
           Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=false;
           Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=false;


           // 18,09,2020 Температура МКЭ задаётся в настройках температуры для amg РУМБА.
            if (Laplas.egddata.iStaticStructural=0) then
          begin
             // Mechanical not active.
             Form_amg_manager.PanelStress1.Visible:=false;
             Form_amg_manager.PanelStress2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogStress.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Stress.Visible:=false;
             Form_amg_manager.CheckBoxStressMatrixPortrait.Visible:=false;

          end;

           Form_amg_manager.ShowModal;
        end;

        if (ComboBoxStaticStructuralSolverSetting.ItemIndex=3) then
        begin  // Алгебраический многосеточный метод amg1r5
           // amg1r5 Ruge and Stueben.
           Formamg1r5Parameters.ShowModal;
        end;

         if (ComboBoxStaticStructuralSolverSetting.ItemIndex=4) then
        begin
           // AMGCL Denis Demidov BiCGStab + samg
           FormAMGCLParameters.ShowModal;
        end;
    end;
    // Графовый метод решения уравнения теплопроводности.
    if (Laplas.egddata.itemper=3) then
    begin
        if (ComboBoxSolverSetting.ItemIndex=5) then
       begin
          // AMGCL Denis Demidov BiCGStab + samg
          FormAMGCLParameters.ShowModal;
       end
       else  if (ComboBoxSolverSetting.ItemIndex=1) then
       begin
          // amg1r5 Ruge and Stueben.
          Formamg1r5Parameters.ShowModal;
       end
       else  if (ComboBoxSolverSetting.ItemIndex=7) then
       begin
          // Вызов настроек РУМБА v.0.14 алгебраического многосеточного метода.

          // cfd not active
           Form_amg_manager.PanelSpeed1.Visible:=false;
           Form_amg_manager.PanelPressure1.Visible:=false;
           Form_amg_manager.PanelSpeed2.Visible:=false;
           Form_amg_manager.PanelPressure2.Visible:=false;
           Form_amg_manager.CheckBoxprintlogSpeed.Visible:=false;
           Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=false;
           Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=false;
           Form_amg_manager.CheckBoxprintlogPressure.Visible:=false;
           Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=false;
           Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=false;

              // Mechanical not active.
             Form_amg_manager.PanelStress1.Visible:=false;
             Form_amg_manager.PanelStress2.Visible:=false;
             Form_amg_manager.CheckBoxprintlogStress.Visible:=false;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Stress.Visible:=false;
             Form_amg_manager.CheckBoxStressMatrixPortrait.Visible:=false;

             // Уравнение теплопередачи решается.
             Form_amg_manager.PanelTemperature1.Visible:=true;
             Form_amg_manager.PanelTemperature2.Visible:=true;
             Form_amg_manager.CheckBoxprintlogTemperature.Visible:=true;
             Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Temperature.Visible:=true;
             Form_amg_manager.CheckBoxTemperatureMatrixPortrait.Visible:=true;

           Form_amg_manager.ShowModal;
       end;
    end;
     if (Laplas.egddata.iStaticStructural=1) then
     begin
       if (ComboBoxSolverSetting.ItemIndex=5) then
       begin
          // AMGCL Denis Demidov BiCGStab + samg
          FormAMGCLParameters.ShowModal;
       end
       else  if (ComboBoxSolverSetting.ItemIndex=1) then
       begin
          // amg1r5 Ruge and Stueben.
          Formamg1r5Parameters.ShowModal;
       end
       else
       begin
          // Вызов настроек РУМБА 0.14
          // cfd not active
          Form_amg_manager.PanelSpeed1.Visible:=false;
          Form_amg_manager.PanelPressure1.Visible:=false;
          Form_amg_manager.PanelSpeed2.Visible:=false;
          Form_amg_manager.PanelPressure2.Visible:=false;
          Form_amg_manager.CheckBoxprintlogSpeed.Visible:=false;
          Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Speed.Visible:=false;
          Form_amg_manager.CheckBoxSpeedMatrixPortrait.Visible:=false;
          Form_amg_manager.CheckBoxprintlogPressure.Visible:=false;
          Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Pressure.Visible:=false;
          Form_amg_manager.CheckBoxPressureMatrixPortrait.Visible:=false;

          // Mechanical is active.
          Form_amg_manager.PanelStress1.Visible:=true;
          Form_amg_manager.PanelStress2.Visible:=true;
          Form_amg_manager.CheckBoxprintlogStress.Visible:=true;
          Form_amg_manager.ComboBoxCFalgorithmandDataStruct_Stress.Visible:=true;
          Form_amg_manager.CheckBoxStressMatrixPortrait.Visible:=true;

          Form_amg_manager.ShowModal;
       end;
    end;
end;

// Управление показом amg mager`a.
procedure TFormSetting.ComboBoxSolverSettingChange(Sender: TObject);
begin


   case ComboBoxSolverSetting.ItemIndex of
     0 : // bicgstab + ilu2
     begin
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=true;
        GroupBox1.Visible:=false; // gmres restart
     end;
     1 :  // amg1r5
     begin
        //  Ruge & Stueben 1987.
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=true;
        GroupBox1.Visible:=true; // gmres restart
     end;
     2 :  // bicgstab + adi
     begin
        if (MeshForm.CheckBoxALICE.Checked) then
        begin
           Laplas.MainMemo.Lines.Add('Do not supported ADI preconditioner in ALICE Mesh.');
           ShowMessage('Do not supported ADI preconditioner in ALICE Mesh.');
           ComboBoxSolverSetting.ItemIndex:=0;
        end;
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
     end;
      3 : // cl_agl_amg_v0_14  РУМБА v0.14
     begin
        // Это более поздняя версия кода в ней
        // ликвидированы positive connections.
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=true;
        GroupBox1.Visible:=true; // gmres restart
     end;
     4 :  // CPU AMGCL Denis Demidov BiCGStab + samg
     begin
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=true; // gmres restart
     end;
     5 :  // GPU AMGCL Denis Demidov BiCGStab + samg
     begin
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=true; // gmres restart
     end;
     6 : // cusp  BiCGStab + samg
     begin
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
     end;
     7 :  // Алгебраический многосеточный метод РУМБА 0.14
     begin
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=true;
        GroupBox1.Visible:=true; // gmres restart
     end;
      8 : // cusp  BiCGStab + samg
     begin
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
     end;
      9 : // ViennaCL 1.7.1  BiCGStab + ILU0
     begin
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
     end;
      10 : //  FGMRES + ILU(lfil)
     begin
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=true;
        GroupBox1.Visible:=true; // gmres restart
     end;
     11 :
     begin   // Chebyshev + amg1r5
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
     12 :
     begin   // Chebyshev + AMGCL CPU
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
     13 :
     begin   // Chebyshev + AMGCL GPU
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
      14 :
     begin   // Chebyshev + Румба v.0.14
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
     15 :
     begin   // Chebyshev + Румба v.0.14
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
     16 :
     begin   // Chebyshev + AMGCL GPU
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
     17 :
     begin   // Chebyshev + AMGCL GPU
       Button_amg_manager.Visible:=true;
       GroupBox_lfil.Visible:=true;
       GroupBox1.Visible:=true; // gmres restart
     end;
   end;
end;

procedure TFormSetting.ComboBoxStaticStructuralSolverSettingChange(
  Sender: TObject);
begin
   if (ComboBoxStaticStructuralSolverSetting.ItemIndex=2) then
   begin  // Алгебраический многосеточный метод РУМБА 0.14
      // Предоставляем доступ к настройкам
      // алгебраического многосеточного метода.
      GroupBox_lfil.Visible:=true;
      GroupBox1.Visible:=true; // gmres restart
      Button_amg_manager.Visible:=true;
   end;
   if (ComboBoxStaticStructuralSolverSetting.ItemIndex=3) then
   begin  // Алгебраический многосеточный метод amg1r5
      // Предоставляем доступ к настройкам
      // алгебраического многосеточного метода.
      GroupBox_lfil.Visible:=true;
      GroupBox1.Visible:=true; // gmres restart
      Button_amg_manager.Visible:=true;
   end;
   if (ComboBoxStaticStructuralSolverSetting.ItemIndex=4) then
   begin
      // AMGCL ddemidov Denis Demidov.
      GroupBox_lfil.Visible:=false;
      GroupBox1.Visible:=true; // gmres restart
      Button_amg_manager.Visible:=true;
   end;
   if (ComboBoxStaticStructuralSolverSetting.ItemIndex=0) then
   begin
      // BiCGStab + ILU(lfil).
      GroupBox_lfil.Visible:=true;
      GroupBox1.Visible:=false; // gmres restart
      Button_amg_manager.Visible:=false;
   end;
   if (ComboBoxStaticStructuralSolverSetting.ItemIndex=1) then
   begin
      // Direct method.
      GroupBox_lfil.Visible:=false;
      GroupBox1.Visible:=false; // gmres restart
      Button_amg_manager.Visible:=false;
   end;
end;

end.
