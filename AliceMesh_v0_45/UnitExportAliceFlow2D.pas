unit UnitExportAliceFlow2D;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormExportAliceFlow2D = class(TForm)
    GroupBox1: TGroupBox;
    Panel1: TPanel;
    CheckBox1: TCheckBox;
    Label1: TLabel;
    Editamplitude: TEdit;
    Label2: TLabel;
    Label3: TLabel;
    Editfrequency: TEdit;
    Label4: TLabel;
    ButtonRun: TButton;
    procedure ButtonRunClick(Sender: TObject);
    procedure CheckBox1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormExportAliceFlow2D: TFormExportAliceFlow2D;

implementation

{$R *.dfm}

procedure TFormExportAliceFlow2D.ButtonRunClick(Sender: TObject);
begin
   Close;
end;

procedure TFormExportAliceFlow2D.CheckBox1Click(Sender: TObject);
begin
   if (CheckBox1.Checked) then
   begin
      Panel1.Visible:=true;
   end
   else
   begin
      Panel1.Visible:=false;
   end;
end;

end.
