unit UnitSquareWave;

// Ќа форму с параметрами желательно добавить красивую картинку чтобы пользователь
// пон€л какие параметры к чему относ€тс€.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls,
  Vcl.Imaging.pngimage;

type
  TFormSquareWave = class(TForm)
    PanelSquareWavegl: TPanel;
    ButtonSquareWaveApply: TButton;
    Labeltau: TLabel;
    Edittau: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Image1: TImage;
    EditQ: TEdit;
    procedure ButtonSquareWaveApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormSquareWave: TFormSquareWave;

implementation

{$R *.dfm}

uses VisualUnit, UnitVariables;

// —читывание параметров термоциклировани€.
procedure TFormSquareWave.ButtonSquareWaveApplyClick(Sender: TObject);
var
   s1 : String;
   k : Integer;
   bOk : Boolean;
   r1 : Real;
begin
  // Laplas.glSTL.iQ:=ComboBoxDutyCycle.ItemIndex+2;
    s1:=Trim(EditQ.Text);
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
   EditQ.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=1.0) then
   begin
      // ‘изический смысл.
      // скважность положительна€ вещественна€ величина больша€ единицы.
      ShowMessage('error user input : Q (Duty cycle) mast be strongly positive > 1.0.');
      bOk:=false;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.iQ:=StrToFloat(EditQ.Text);
      EditQ.Color:=clWhite;
   end
   else
   begin
      EditQ.Text:=FloatToStr(Laplas.glSTL.iQ);
      EditQ.Color:=clRed;
   end;

   // «десь нужно прописать защиту от некорректного ввода.
   s1:=Trim(Edittau.Text);
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
   Edittau.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=0.0) then
   begin
      // ‘изический смысл.
      // длительность импульса положительна€ вещественна€ величина.
      ShowMessage('error user input : tau mast be strongly positive.');
      bOk:=false;
   end;
   if (bOk) then
   begin
      Laplas.glSTL.tau:=StrToFloat(Edittau.Text);
      Edittau.Color:=clWhite;
   end
   else
   begin
      Edittau.Text:=FloatToStr(Laplas.glSTL.tau);
      Edittau.Color:=clRed;
   end;
   if bOk then Close;
end;

end.
