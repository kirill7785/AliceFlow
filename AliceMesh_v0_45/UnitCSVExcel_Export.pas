unit UnitCSVExcel_Export;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormCSV_Excel_export = class(TForm)
    RadioGroup1: TRadioGroup;
    ButtonExport: TButton;
    procedure ButtonExportClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormCSV_Excel_export: TFormCSV_Excel_export;

implementation

{$R *.dfm}

procedure TFormCSV_Excel_export.ButtonExportClick(Sender: TObject);
begin
   Close;
end;

end.
