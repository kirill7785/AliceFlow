unit UnitSetHpolyinReadOBJFiles;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormHpolyinOBJFiles = class(TForm)
    ButtonApply: TButton;
    Label1: TLabel;
    Labelreadblockname: TLabel;
    LabeliPlane: TLabel;
    LabelreadiPlaneObj2: TLabel;
    Label2: TLabel;
    EditHpoly: TEdit;
    Label3: TLabel;
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormHpolyinOBJFiles: TFormHpolyinOBJFiles;

implementation

{$R *.dfm}

procedure TFormHpolyinOBJFiles.ButtonApplyClick(Sender: TObject);
var
  iscan : Integer;
  sub : String;
begin
   sub:= EditHpoly.Text;
   for iscan:=1 to length(sub) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (sub[iscan]=',') then
         begin
            sub[iscan]:='.';
         end;
      end;
   end;

   EditHpoly.Text:=Trim(sub);

   Close;
end;

end.
