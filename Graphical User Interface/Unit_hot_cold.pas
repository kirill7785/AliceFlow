unit Unit_hot_cold;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.Imaging.pngimage,
  Vcl.ExtCtrls;

type
  TForm_hot_cold = class(TForm)
    Image_double_linear: TImage;
    Button1: TButton;
    Label1: TLabel;
    Edit_on_time: TEdit;
    Label_time_union: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form_hot_cold: TForm_hot_cold;

implementation

{$R *.dfm}

uses UnitVariables, VisualUnit;

procedure TForm_hot_cold.Button1Click(Sender: TObject);
var
   s1 : String;
   k : Integer;
   bOk : Boolean;
   r1 : Real;
begin
    s1:=Trim(Edit_on_time.Text);
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
   Edit_on_time.Text:=Trim(s1);
   bOk:=true;
   if bOk then r1:=FormVariables.my_real_convert(s1,bOk);
   if (r1<=0.0) then
   begin
      // Физический смысл.
      // скважность положительная вещественная величина большая единицы.
      ShowMessage('error user input :  mast be strongly positive > 1.0.');
      bOk:=false;
   end
   else
   begin
       Laplas.glSTL.on_time_double_linear:=r1;
   end;
   if (bOk) then
   begin
      Close;
   end;
end;

end.
