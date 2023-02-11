unit UnitMoveSTL_position;
//  оординаты смещени€ позиций всех STL-геометрических объектов.
// 29.11.2022

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormMoveSTLposition = class(TForm)
    Label1: TLabel;
    Label2: TLabel;
    EditX: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    EditY: TEdit;
    Label5: TLabel;
    Label6: TLabel;
    EditZ: TEdit;
    Label7: TLabel;
    Button1: TButton;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormMoveSTLposition: TFormMoveSTLposition;

implementation

{$R *.dfm}

uses VisualUnit, UnitVariables;

procedure TFormMoveSTLposition.Button1Click(Sender: TObject);
var
   bOk : Boolean;
   s1, s2, s3  : String;
   r1, r2, r3 : Real;
begin
   // инициализаци€ :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;

    bOk:=true; // признак правильности ввода

    s1:=EditX.Text;
    s2:=EditY.Text;
    s3:=EditZ.Text;

    if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
    if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
    if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учЄтом подстановки значений переменных.


      if (bOk) then
      begin
         //Laplas.move_STL_position.x:=StrToFloat(EditX.Text);
         //Laplas.move_STL_position.y:=StrToFloat(EditY.Text);
         //Laplas.move_STL_position.z:=StrToFloat(EditZ.Text);

         Laplas.move_STL_position.x:=r1;
         Laplas.move_STL_position.y:=r2;
         Laplas.move_STL_position.z:=r3;

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
