unit UnitRenameVariable;
// Переименовывает пользовательскую переменную в проекте.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormRenameVar = class(TForm)
    GroupBox1: TGroupBox;
    ComboBox1: TComboBox;
    LabelRenamecandidate: TLabel;
    Label1: TLabel;
    EditNewName: TEdit;
    Button1: TButton;
    procedure FormCreate(Sender: TObject);
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    bOk_rename : Boolean;
  end;

var
  FormRenameVar: TFormRenameVar;

implementation

{$R *.dfm}

uses VisualUnit;

procedure TFormRenameVar.Button1Click(Sender: TObject);
var
   bfound_variable : Boolean;
   j_var : Integer;
   sub : String;
begin
   EditNewName.Color:=clwhite;
   if (length(Trim(EditNewName.Text))>1) then
   begin
       if (Pos('$',Trim(EditNewName.Text))>0) then
       begin
           // Новое имя не совпадает с предыдущими именами переменных.
            bfound_variable:=false;
            sub:= Trim(EditNewName.Text);
            for j_var := 0 to Laplas.ivar-1 do
            begin
               if ((length(Trim(sub))=length(Trim(Laplas.parametric[j_var].svar)))and(pos(Trim(sub),Trim(Laplas.parametric[j_var].svar))=1))  then
               begin
                  // Найдено полное совпадение.
                  bfound_variable:=true;
               end;
            end;
           if (not(bfound_variable)) then
           begin
              bOk_rename:=true;
              Close;
           end;
       end;
   end;
   if (not(bOk_rename)) then
   begin
      EditNewName.Color:=clred;
   end;
end;

procedure TFormRenameVar.FormCreate(Sender: TObject);
begin
   bOk_rename:=false;
end;

end.
