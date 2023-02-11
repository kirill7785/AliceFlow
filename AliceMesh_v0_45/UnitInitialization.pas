unit UnitInitialization;
// Задает начальную скорость. Задает скорость
// при решении уравнения конвекции-диффузии.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormSpeedInitialization = class(TForm)
    PanelInitializationMain: TPanel;
    EditVx: TEdit;
    EditVy: TEdit;
    EditVz: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Button1: TButton;
    procedure Button1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormSpeedInitialization: TFormSpeedInitialization;

implementation

{$R *.dfm}

procedure TFormSpeedInitialization.Button1Click(Sender: TObject);
var
   s : String;
begin
    if (FormatSettings.DecimalSeparator=',') then
    begin
        // заменить все точки в FormSpeedInitialization на запятые.
        s:=FormSpeedInitialization.EditVx.Text;
        FormSpeedInitialization.EditVx.Text:=StringReplace(s,'.',',',[rfReplaceAll]);

        s:=FormSpeedInitialization.EditVy.Text;
        FormSpeedInitialization.EditVy.Text:=StringReplace(s,'.',',',[rfReplaceAll]);

        s:=FormSpeedInitialization.EditVz.Text;
        FormSpeedInitialization.EditVz.Text:=StringReplace(s,'.',',',[rfReplaceAll]);

    end;
     if (FormatSettings.DecimalSeparator='.') then
    begin
        // заменить все точки в FormSpeedInitialization на запятые.
        s:=FormSpeedInitialization.EditVx.Text;
        FormSpeedInitialization.EditVx.Text:=StringReplace(s,',','.',[rfReplaceAll]);

        s:=FormSpeedInitialization.EditVy.Text;
        FormSpeedInitialization.EditVy.Text:=StringReplace(s,',','.',[rfReplaceAll]);

        s:=FormSpeedInitialization.EditVz.Text;
        FormSpeedInitialization.EditVz.Text:=StringReplace(s,',','.',[rfReplaceAll]);

    end;

    // Нужна проверка на корректность преобразования в вещественное число из строкового типа.
     Close();
end;

end.
