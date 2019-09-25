unit UnitSquareWaveAPPARAT;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls,
  Vcl.Imaging.pngimage, Vcl.AppEvnts;

type
  TFormAPPARAT_Square_Wave = class(TForm)
    Panel1: TPanel;
    Label1: TLabel;
    Editmultiplyer: TEdit;
    Labeltau1: TLabel;
    Edittau1: TEdit;
    Label2: TLabel;
    Labeltau2: TLabel;
    Edittau2: TEdit;
    Label3: TLabel;
    ButtonApplySquareWaveApparat: TButton;
    Labeltaupause: TLabel;
    Edit_tau_pause: TEdit;
    Label4: TLabel;
    Label5: TLabel;
    ComboBoxnumbercycle: TComboBox;
    Label6: TLabel;
    EditPeriod: TEdit;
    Label7: TLabel;
    Image1: TImage;
    ApplicationEvents1: TApplicationEvents;
    procedure ButtonApplySquareWaveApparatClick(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormAPPARAT_Square_Wave: TFormAPPARAT_Square_Wave;

implementation

{$R *.dfm}

uses
    VisualUnit, UnitVariables;

// «апрет форме сворачиватьс€.
procedure TFormAPPARAT_Square_Wave.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

procedure TFormAPPARAT_Square_Wave.ButtonApplySquareWaveApparatClick(Sender: TObject);
var
   s1 : String;
   k : Integer;
   bOk : Boolean;
   r1 : Real;
begin
   s1:=Trim(Editmultiplyer.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   Editmultiplyer.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
    if ((r1<=0.0)or(r1>=1.0)) then
   begin
      // ‘изический смысл.
      // длительность импульса положительна€ вещественна€ величина.
      ShowMessage('error user input : multiplyer mast be in open interval (0.0..1.0).');
      bOk:=false;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.m1:=StrToFloat(Editmultiplyer.Text);
      Editmultiplyer.Color:=clWhite;
   end
   else
   begin
      Editmultiplyer.Text:=FloatToStr(Laplas.glSTL.m1);
      Editmultiplyer.Color:=clRed;
   end;

   s1:=Trim(Edittau1.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   Edittau1.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=0.0) then
   begin
      // ‘изический смысл.
      // длительность импульса положительна€ вещественна€ величина.
      ShowMessage('error user input : tau1 mast be strongly positive.');
      bOk:=false;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.tau1:=StrToFloat(Edittau1.Text);
      Edittau1.Color:=clWhite;
   end
   else
   begin
      Edittau1.Text:=FloatToStr(Laplas.glSTL.tau1);
      Edittau1.Color:=clRed;
   end;

    s1:=Trim(Edittau2.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   Edittau2.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=0.0) then
   begin
      // ‘изический смысл.
      // длительность импульса положительна€ вещественна€ величина.
      ShowMessage('error user input : tau2 mast be strongly positive.');
      bOk:=false;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.tau2:=StrToFloat(Edittau2.Text);
      Edittau1.Color:=clWhite;
   end
   else
   begin
      Edittau2.Text:=FloatToStr(Laplas.glSTL.tau2);
      Edittau2.Color:=clRed;
   end;

    s1:=Trim(Edit_tau_pause.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   Edit_tau_pause.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=0.0) then
   begin
      // ‘изический смысл.
      // длительность импульса положительна€ вещественна€ величина.
      ShowMessage('error user input : tau_pause mast be strongly positive.');
      bOk:=false;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.tau_pause:=StrToFloat(Edit_tau_pause.Text);
      Edittau1.Color:=clWhite;
   end
   else
   begin
      Edit_tau_pause.Text:=FloatToStr(Laplas.glSTL.tau_pause);
      Edit_tau_pause.Color:=clRed;
   end;

     s1:=Trim(EditPeriod.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditPeriod.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=0.0) then
   begin
      // ‘изический смысл.
      // длительность импульса положительна€ вещественна€ величина.
      ShowMessage('error user input : Period mast be strongly positive.');
      bOk:=false;
   end;
   if (bOk) then
   begin
      if (r1<=Laplas.glSTL.n*(2*Laplas.glSTL.tau1+Laplas.glSTL.tau2+Laplas.glSTL.tau_pause)) then
      begin
         bOk:=false;
         ShowMessage('error user input : Period mast be > number_cycle*(2*tau1+tau2+tau_pause).');
      end;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.T:=StrToFloat(EditPeriod.Text);
      EditPeriod.Color:=clWhite;
   end
   else
   begin
      EditPeriod.Text:=FloatToStr(Laplas.glSTL.T);
      EditPeriod.Color:=clRed;
   end;

   Laplas.glSTL.n:=ComboBoxnumbercycle.ItemIndex+1;

   if (bOk) then
   begin
      Close;
   end;
end;

end.
