unit AddVariableUnit;
// Добавляет переменную в проект,
// Если пользователь забыл её ввести.


interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TAddVariableForm = class(TForm)
    pnlbase: TPanel;
    lblnamevar: TLabel;
    lblname: TLabel;
    lblvalue: TLabel;
    edtvalue: TEdit;
    btnApply: TButton;
    procedure btnApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  AddVariableForm: TAddVariableForm;

implementation

uses UnitVariables, VisualUnit;

{$R *.dfm}

// По нажатию на кнопку Apply мы должны ввести новую переменную в
// список уже существующих переменных.
procedure TAddVariableForm.btnApplyClick(Sender: TObject);
var
  // переменные которые понадобяться для добавления новой переменной.
    i, code : Integer;
    b : Boolean;
    s : string;
    c : Real;
    parametric_copy : array of TmyVariable;
begin

     // Вводит переменные в проект.
     // переменные должны идти последовательно сверху вниз.
     // переменные не удовлетворяющие этому критерию исключаются из рассмотрения.
     // имя переменной обязательно должно начинаться со знака $.
     b:=true;
     i:=0;
     while (b and (i<FormVariables.StringGridVariables.RowCount-2)) do
     begin
        s:=FormVariables.StringGridVariables.Cells[1,i+1];
        if ((length(s)>=2) and (s[1]='$')) then b:=true
         else b:=false;
        if (b) then i:=i+1; // переходим к следующей переменной
     end;
     Laplas.ivar:=i; // количество переменных.
     SetLength(Laplas.parametric,Laplas.ivar);
     for i:=0 to Laplas.ivar-1 do
     begin
        Laplas.parametric[i].svar:=Trim(FormVariables.StringGridVariables.Cells[1,i+1]);
        Laplas.parametric[i].sval:=Trim(FormVariables.StringGridVariables.Cells[2,i+1]);
     end;
     // Все переменные хранятся в таблице  StringGridVariables

     // Все последующине переменные игнорируются и их надо стереть
     for i:=Laplas.ivar+1 to FormVariables.StringGridVariables.RowCount-1 do
     begin
        FormVariables.StringGridVariables.Cells[1,i]:='';
        FormVariables.StringGridVariables.Cells[2,i]:='';
     end;

     lblname.Caption:=Trim(lblname.Caption);
     if (lblname.Caption[1]<>'$') then
     begin
        lblname.Caption:='$'+lblname.Caption; // добавление доллара как символа переменной.
     end;
     // Добавление новой переменной.
     val(edtvalue.Text,c,code); // пытаемся получить число
     if code=0 then // константа
     begin
        SetLength(parametric_copy,Laplas.ivar);
        for i:=0 to Laplas.ivar-1 do
        begin
           parametric_copy[i]:=Laplas.parametric[i];
        end;
        Inc(Laplas.ivar); // увеличиваем количество переменных на единицу.
        FormVariables.StringGridVariables.Cells[1,Laplas.ivar]:=Trim(lblname.Caption);
        FormVariables.StringGridVariables.Cells[2,Laplas.ivar]:=Trim(edtvalue.Text);
        SetLength(Laplas.parametric,Laplas.ivar);
        for i:=0 to Laplas.ivar-2 do
        begin
           Laplas.parametric[i]:=parametric_copy[i];
        end;
        Laplas.parametric[Laplas.ivar-1].svar:=Trim(lblname.Caption);
        Laplas.parametric[Laplas.ivar-1].sval:=Trim(edtvalue.Text);
        Close; // уходим из данной формы.
     end
      else
     begin
        // откат назад и всё по новой.
        edtvalue.Text:='0.0';
        FormVariables.StringGridVariables.Cells[1,Laplas.ivar]:='';
        FormVariables.StringGridVariables.Cells[2,Laplas.ivar]:='';
        dec(Laplas.ivar);
        ShowMessage('Error! Please, enter value variable correctly...');
     end;                                                               
end;

end.
