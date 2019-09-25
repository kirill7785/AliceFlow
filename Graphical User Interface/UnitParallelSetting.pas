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
  Unitamg1r5Parameters;


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
       Form_amg_manager.ShowModal;
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
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
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
     4 :  // cusp   BiCGStab+AINV (NS Bridson)
     begin
        Button_amg_manager.Visible:=false;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
     end;
     5 :  // AMGCL Denis Demidov BiCGStab + samg
     begin
        Button_amg_manager.Visible:=true;
        GroupBox_lfil.Visible:=false;
        GroupBox1.Visible:=false; // gmres restart
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
     begin   // CUSP BiCGStab + AINV
       Button_amg_manager.Visible:=false;
       GroupBox_lfil.Visible:=false;
       GroupBox1.Visible:=false; // gmres restart
     end;
   end;
end;

procedure TFormSetting.ComboBoxStaticStructuralSolverSettingChange(
  Sender: TObject);
begin
   if (ComboBoxStaticStructuralSolverSetting.ItemIndex=2) then
   begin
      // Предоставляем доступ к настройкам алгебраического многосеточного метода.
      Button_amg_manager.Visible:=true;
   end;
end;

end.
