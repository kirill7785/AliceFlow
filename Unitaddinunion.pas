unit Unitaddinunion;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, Vcl.AppEvnts;

type
  TForm_formirateunion = class(TForm)
    grp1: TGroupBox;
    cbbselectnameunion: TComboBox;
    btnApply: TButton;
    ApplicationEvents1: TApplicationEvents;
    procedure btnApplyClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
  private
    { Private declarations }
  public
    { Public declarations }
    iselectedunion : Integer; // запоминаем выбранный в последний раз union.
  end;

var
  Form_formirateunion: TForm_formirateunion;

implementation

{$R *.dfm}

 // Запрет форме сворачиваться.
procedure TForm_formirateunion.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
     if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

procedure TForm_formirateunion.btnApplyClick(Sender: TObject);
begin
   // вся обработка будет в модуле VisualUnit.
   Close;
end;

procedure TForm_formirateunion.FormCreate(Sender: TObject);
begin
   iselectedunion:=0;
end;

end.
