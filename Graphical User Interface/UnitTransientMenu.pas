unit UnitTransientMenu;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormUnsteady = class(TForm)
    RadioGroup1: TRadioGroup;
    PanelTime: TPanel;
    LabelEndTime: TLabel;
    EditTime: TEdit;
    LabelTimeUnion: TLabel;
    ComboBoxTimeStep: TComboBox;
    ButtonTimeStepLaw: TButton;
    GroupBox1: TGroupBox;
    CheckBoxdonttec360: TCheckBox;
    CheckBoxonlysolidvisible: TCheckBox;
    CheckBoxreconstruct: TCheckBox;
    CheckBoxAnimationFields: TCheckBox;
    CheckBoxCylinderToPrism: TCheckBox;
    ButtonCalc: TButton;
    GroupBoxNumberIterationsSimpleAlgorithm: TGroupBox;
    ComboBoxNumberIterationsSIMPLEAlgoritm: TComboBox;
    procedure ButtonCalcClick(Sender: TObject);
    procedure RadioGroup1Click(Sender: TObject);
    procedure ButtonTimeStepLawClick(Sender: TObject);
    procedure ComboBoxTimeStepChange(Sender: TObject);
    procedure FormClose(Sender: TObject; var Action: TCloseAction);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    bRunOk : Boolean;
  end;

var
  FormUnsteady: TFormUnsteady;

implementation

{$R *.dfm}

uses UnitVariables, UnitSquareWave, VisualUnit, UnitSquareWaveAPPARAT,
  Unit_hot_cold;

procedure TFormUnsteady.ButtonCalcClick(Sender: TObject);
var
   s : String;
   bOk : Boolean;
begin
    if (length(EditTime.Text)>0) then
    begin
        s:=EditTime.Text;
        bOk:=true;
        FormVariables.my_real_convert(s,bOk);
        if bOk then
        begin
           EditTime.Color:=clWhite;
           Close;
        end
        else
        begin
           EditTime.Color:=clRed;
        end;
    end
    else
    begin
       EditTime.Color:=clRed;
    end;
    bRunOk:=true;
end;


// –едактирвание параметров импульсного режима.
procedure TFormUnsteady.ButtonTimeStepLawClick(Sender: TObject);
begin
   if (ComboBoxTimeStep.ItemIndex=1) then
   begin
     // Square Wave Law.
      FormSquareWave.EditQ.Text:=FloatToStr(Laplas.glSTL.iQ);
      FormSquareWave.Edittau.Text:=FloatToStr(Laplas.glSTL.tau);
      FormSquareWave.ShowModal;
   end;
   if (ComboBoxTimeStep.ItemIndex=2) then
   begin
      // Square Wave APPARAT.
      FormAPPARAT_Square_Wave.EditPeriod.Text:=FloatToStr(Laplas.glSTL.T);
      FormAPPARAT_Square_Wave.Editmultiplyer.Text:=FloatToStr(Laplas.glSTL.m1);
      FormAPPARAT_Square_Wave.Edittau1.Text:=FloatToStr(Laplas.glSTL.tau1);
      FormAPPARAT_Square_Wave.Edittau2.Text:=FloatToStr(Laplas.glSTL.tau2);
      FormAPPARAT_Square_Wave.Edit_tau_pause.Text:=FloatToStr(Laplas.glSTL.tau_pause);
      FormAPPARAT_Square_Wave.ComboBoxnumbercycle.ItemIndex:=Laplas.glSTL.n-1;
      FormAPPARAT_Square_Wave.ShowModal;
   end;
   // double linear time step
   if (ComboBoxTimeStep.ItemIndex=3) then
   begin
      Form_hot_cold.Edit_on_time.Text:=FloatToStr(Laplas.glSTL.on_time_double_linear);
      Form_hot_cold.ShowModal;
   end;
end;

procedure TFormUnsteady.ComboBoxTimeStepChange(Sender: TObject);
begin
   case ComboBoxTimeStep.ItemIndex of
     0 : begin
       ButtonTimeStepLaw.Visible:=false;
     end;
     1 : begin
         ButtonTimeStepLaw.Visible:=true;
     end;
     2 : begin
        ButtonTimeStepLaw.Visible:=true;
     end;
     3 : begin
        ButtonTimeStepLaw.Visible:=true;
     end;
   end;
end;

procedure TFormUnsteady.FormClose(Sender: TObject; var Action: TCloseAction);
begin
  if (not(bRunOk)) then
  begin
     bRunOk:=false;
  end;
end;

procedure TFormUnsteady.FormCreate(Sender: TObject);
begin
   bRunOk:=false;
end;

// ѕоказывать ли задание временного интервала дл€
// нестационарного моделировани€.
procedure TFormUnsteady.RadioGroup1Click(Sender: TObject);
begin
    case RadioGroup1.ItemIndex of
       0 : begin
              PanelTime.Visible:=false;
           end;
       1 : begin
              // 31.08.2019
              if (Laplas.egddata.myflmod[0].iflow=0) then
              begin
                 // „иста€ теплопередача
                 PanelTime.Visible:=True;
              end
               else
              begin
                 // √идродинамика не может быть нестационарной.
                 RadioGroup1.ItemIndex:=0;
                 Laplas.MainMemo.Lines.Add('WARNING!!! unsteady cfd calculation not support.');
                 Application.MessageBox('unsteady cfd calculation not support.','WARNING',MB_OK)
              end;
           end;
    end;
end;

end.
