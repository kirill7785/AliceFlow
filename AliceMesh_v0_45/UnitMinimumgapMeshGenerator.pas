unit UnitMinimumgapMeshGenerator;
// Настройка сеточного генератора минимальный зазор между двумя сеточными линиями, всё что меньше будет игнорироваться.
// 29.11.2022

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormMinimum_gap = class(TForm)
    Label1: TLabel;
    Label2: TLabel;
    EditX: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Button1: TButton;
    EditY: TEdit;
    EditZ: TEdit;
    Label6: TLabel;
    Label7: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormMinimum_gap: TFormMinimum_gap;

implementation

{$R *.dfm}

uses VisualUnit, UnitVariables;

procedure TFormMinimum_gap.Button1Click(Sender: TObject);
var
   bOk : Boolean;
   s1, s2, s3  : String;
   r1, r2, r3 : Real;
begin

   // инициализация :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;

    bOk:=true; // признак правильности ввода

    s1:=EditX.Text;
    s2:=EditY.Text;
    s3:=EditZ.Text;

    if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
    if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
    if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учётом подстановки значений переменных.


      if (bOk) then
      begin

         //Laplas.minimum_gap.x:=StrToFloat(EditX.Text);
         //Laplas.minimum_gap.y:=StrToFloat(EditY.Text);
         //Laplas.minimum_gap.z:=StrToFloat(EditZ.Text);

         Laplas.minimum_gap.x:=r1;
         Laplas.minimum_gap.y:=r2;
         Laplas.minimum_gap.z:=r3;

         EditX.Color:=clwhite;
         EditY.Color:=clwhite;
         EditZ.Color:=clwhite;

      end
      else
      begin
         EditX.Color:=clred;
         EditY.Color:=clred;
         EditZ.Color:=clred;
      end;
end;

end.
